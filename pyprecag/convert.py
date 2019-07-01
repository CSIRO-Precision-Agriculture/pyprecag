# -*- coding: utf-8 -*-
"""
Format conversion routines.
"""
from  datetime import timedelta
import inspect
import logging
import math
import os
import time

import numpy as np

import pandas as pd
from pandas.core.dtypes.common import is_string_dtype

from osgeo import ogr, gdal
import geopandas
from geopandas import GeoDataFrame, GeoSeries

import rasterio
from rasterio import features

from fiona import collection as fionacoll
from fiona.crs import from_epsg
from shapely.geometry import Point, mapping, shape, LineString

from . import crs as pyprecag_crs
from . import TEMPDIR, config
from .describe import CsvDescribe, predictCoordinateColumnNames, VectorDescribe, \
    save_geopandas_tofile

from .errors import GeometryError
from .raster_ops import raster_snap_extent, create_raster_transform

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def numeric_pixelsize_to_string(pixelsize_metres):
    """
    Convert a numeric pixel size to a string with the appropriate units for use in file names.
    For example 0.40 is 40cm,  0.125 is 125mm, 1.5km is 1500m,  2.0 is 2m, 6000 is 6km
    Args:
        pixelsize_metres (numeric): The pixel size in meters.

    Returns (str): a string representation of the pixel size
    """
    if isinstance(pixelsize_metres, str):
        raise ValueError('Pixel size should be numeric (int, float etc)')

    if pixelsize_metres % 1000.0 == 0:  # or km
        result = '{}km'.format(int(pixelsize_metres / 1000))
    elif pixelsize_metres % 1 == 0:  # is in m
        result = '{}m'.format(int(pixelsize_metres))
    elif (pixelsize_metres * 100.0) % 1 == 0:  # is in cm
        result = '{}cm'.format(int(pixelsize_metres * 100))
    else:  # or mm
        result = '{}mm'.format(int(pixelsize_metres * 1000))
    return result


def point_to_point_bearing(pt1, pt2):
    """Calculate bearing from north between two points """
    # variation of  https://community.esri.com/thread/27393

    if isinstance(pt1, tuple):
        pt1 = Point(pt1)

    if isinstance(pt2, tuple):
        pt2 = Point(pt2)

    bearing = 180 + math.atan2((pt2.x - pt1.x), (pt2.y - pt1.y)) * (180 / math.pi)
    return bearing


def line_bearing(line):
    """ The bearing on the line between the start and end nodes/points"""
    if not isinstance(line, LineString):
        line = line.geometry
    return point_to_point_bearing(line.coords[-1], line.coords[0])


def deg_to_16_compass_pts(degrees):
    """
    Convert a compass bearing from north using 16 cardinal directions
    Args:
        degrees (float): Degrees from north
    Returns:
        direction(str) : letters representing the compass point
    source: https://stackoverflow.com/a/7490772
    """

    # fix if a negative number is used ie -45 turns into 315
    if degrees < 0:
        degrees = 360.0 + degrees
    val = int((degrees / 22.5) + .5)

    arr = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW",
           "NW", "NNW"]
    return arr[(val % 16)]


def deg_to_8_compass_pts(degrees):
    """
    Convert a compass bearing from north using 8 cardinal directions
    Args:
        degrees (float): Degrees from north
    Returns:
        direction(str) : letters representing the compass point

    source: variation of https://stackoverflow.com/a/7490772
    """
    # fix if a  negative number is used
    if degrees < 0:
        degrees = 360.0 + degrees

    val = int((degrees / 45.0) + .5)

    arr = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]

    return arr[(val % 8)]


def drop_z(geom):
    if isinstance(geom, LineString):
        return LineString([xy[0:2] for xy in list(geom.coords)])
    if isinstance(geom, Point):
        return Point([geom.coords[0:2]])


def text_rotation(pt1, pt2):
    # calculate an angle to use for labelling in matplotlib plots
    # https://stackoverflow.com/a/35825527
    rotation = np.rad2deg(np.arctan2(pt2.y - pt1.y, pt2.x - pt1.x))
    # need to adjust slighlt to avoid upside down labels
    if rotation < -100:
        rotation = np.rad2deg(np.arctan2(pt1.y - pt2.y, pt1.x - pt2.x))
    elif rotation >= 170:
        rotation = rotation - 180
    return rotation


def convert_polygon_to_grid(in_shapefilename,
                            out_rasterfilename,
                            pixel_size,
                            nodata_val=-9999,
                            snap_extent_to_pixel=True,
                            overwrite=True):
    """ Convert a polygon shapefile to a raster for a set pixel size.

    source : http://ceholden.github.io/open-geo-tutorial/python/chapter_4_vector.html

    Args:
        in_shapefilename (str): A polygon shapefile to rasterise
        out_rasterfilename (str): the output raster name
        pixel_size (float): the pixel size for the raster
        nodata_val (int):the value to use for nodata
        snap_extent_to_pixel (bool):round the extent coordinates to be divisible by the pixel size
        overwrite (bool): if true overwrite existing output file

    Returns:
        None:
    """

    start_time = time.time()

    if not os.path.exists(in_shapefilename):
        raise IOError("Invalid path: {}".format(in_shapefilename))

    shp_ds = ogr.Open(in_shapefilename)
    shp_layer = shp_ds.GetLayer()

    x_min, x_max, y_min, y_max = shp_layer.GetExtent()

    LOGGER.info('\tVector Extent: LL:{} {}'.format(x_min, y_min))
    LOGGER.info('\t               UR:{} {}'.format(x_max, y_max))

    # Snap output grids to a multiple of the grid size, allowing adjacent blocks to align nicely.
    if snap_extent_to_pixel:
        x_min, y_min, x_max, y_max = raster_snap_extent(x_min, y_min, x_max, y_max, pixel_size)
        LOGGER.info('\tRaster Extent: Snap is {}'.format(snap_extent_to_pixel))
        LOGGER.info('\t           LL:{} {}'.format(x_min, y_min))
        LOGGER.info('\t           UR:{} {}'.format(x_max, y_max))

    x_cols = int((x_max - x_min) / pixel_size)
    y_rows = int((y_max - y_min) / pixel_size)
    LOGGER.info('\t           Cols: {}'.format(x_cols))
    LOGGER.info('\t           Rows: {}'.format(y_rows))

    # for pixel data types see http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4
    tif_drv = gdal.GetDriverByName('GTiff')
    if tif_drv is None:
        # If the above Fails then here is an alternate version.
        tif_drv = Find_GDALDriverByName('GTiff')

    raster_ds = tif_drv.Create(out_rasterfilename, x_cols, y_rows, 1, gdal.GDT_Int16)

    # Define spatial reference
    raster_ds.SetProjection(shp_layer.GetSpatialRef().ExportToWkt())
    ras_y_max = y_min + (y_rows * pixel_size)
    raster_ds.SetGeoTransform((x_min, pixel_size, 0, ras_y_max, 0, -pixel_size))

    r_band = raster_ds.GetRasterBand(1)
    r_band.SetNoDataValue(nodata_val)
    r_band.Fill(nodata_val)

    # Rasterize the shapefile layer to our new dataset
    # options = ['ALL_TOUCHED=TRUE','ATTRIBUTE=id']
    err = gdal.RasterizeLayer(raster_ds,  # output to our new dataset
                              [1],  # output to our new datasets first band
                              shp_layer,  # rasterize this layer
                              None, None,  # ignore transformations since we're in same projection
                              burn_values=[1],  # burn value 1
                              options=['ALL_TOUCHED=FALSE'])

    if err != 0:
        LOGGER.exception("Error rasterizing layer: {}".format(err), True)

    else:
        LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                                   dur=timedelta(seconds=time.time() - start_time)))

    # close the data sources
    r_band = None
    raster_ds = None
    shp_layer = None
    shp_ds = None

    return


def convert_grid_to_vesper(in_rasterfilename, out_vesperfilename):
    """ Convert a Grid the format compatible with vesper.

       The vesper format contains the xy coordinates for the center of each pixel with one space
       before the x and three spaces between the x and y values.

       The coordinates are written as rows from the Upper Left corner of the raster.

       Args:
           in_rasterfilename (str):An input raster file
           out_vesperfilename (str): Save the output vesper file to this file.

    """
    start_time = time.time()

    if not os.path.exists(in_rasterfilename):
        raise IOError("Invalid path: {}".format(in_rasterfilename))

    import rasterio  # for some reason this works better here than at the top
    with rasterio.open(os.path.normpath(in_rasterfilename)) as oRasterDS:
        band1 = oRasterDS.read(1)
        with open(out_vesperfilename, 'w') as vesFile:
            for iRowY in range(oRasterDS.height):  # get starting row
                for iColX in (range(oRasterDS.width)):  # process through cols
                    val = band1[iRowY, iColX]
                    if val > oRasterDS.nodatavals:
                        # Get XY coordinates for current row/col
                        # takes rows,cols and returns xy coord for cell center
                        xy = oRasterDS.xy(iRowY, iColX)
                        vesFile.write(' {}   {}\n'.format(xy[0], xy[1]))

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=timedelta(seconds=time.time() - start_time)))
    return


def add_point_geometry_to_dataframe(in_dataframe, coord_columns=None,
                                    coord_columns_epsg=0, out_epsg=0):
    """ Add Point geometry to a pandas or geopandas data frame through sets of coordinate columns.
       If required, it will also reproject to the chosen EPSG coordinate system

    Args:
        in_dataframe (DataFrame): The input dataframe. This could be a pandas or geopandas dataframe
        coord_columns (list) : A list of the columns representing coordinates. If None, it will be
                      predicted.
        coord_columns_epsg (int):  The EPSG number representing the coordinates inside the csv file
        out_epsg (int): The EPSG number representing the required output coordinates system.
                        Use  0 for no projection
    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.

    Modified from: https://gis.stackexchange.com/a/182234

    """
    start_time = time.time()

    if not isinstance(in_dataframe, (geopandas.GeoDataFrame, pd.DataFrame)):
        raise TypeError("Input points should be a geopandas dataframe")

    for argCheck in [('coord_columns_epsg', coord_columns_epsg), ('out_epsg', out_epsg)]:
        if not isinstance(argCheck[1], (int, long)):
            raise TypeError('{} must be a integer.'.format(argCheck[0]))

    if coord_columns is None:
        coord_columns = []

    if not isinstance(coord_columns, list):
        raise TypeError('Coordinate columns should be a list.'.format(coord_columns))

    elif len(coord_columns) > 0:
        for eaFld in coord_columns:
            if eaFld not in in_dataframe.columns:
                raise TypeError('Column {} does not exist'.format(eaFld))
        x_column, y_column = coord_columns

    else:
        x_column, y_column = predictCoordinateColumnNames(in_dataframe.columns)

    if 'FID' not in in_dataframe.columns:
        # insert feature id as column as first column and populate it using row number.
        in_dataframe.insert(0, 'FID', in_dataframe.index)

    # combine lat and lon column to a shapely Point() object
    in_dataframe['geometry'] = in_dataframe.apply(
        lambda x: Point((float(x[x_column]), float(x[y_column]))), axis=1)

    # Now, convert the pandas DataFrame into a GeoDataFrame. The geopandas constructor expects a
    # geometry column which can consist of shapely geometry objects, so the column we created is
    # just fine:
    gdf_csv = geopandas.GeoDataFrame(in_dataframe, geometry='geometry',
                                     crs=from_epsg(coord_columns_epsg))

    # drop the original geometry columns to avoid confusion
    gdf_csv.drop([x_column, y_column], axis=1, inplace=True)

    gdf_crs = pyprecag_crs.crs()
    gdf_crs.getFromEPSG(coord_columns_epsg)

    if out_epsg == -1:
        xmin, ymin, _, _ = gdf_csv.total_bounds
        out_epsg = pyprecag_crs.getProjectedCRSForXY(xmin, ymin, coord_columns_epsg).epsg_number

    if out_epsg > 0:
        gdf_csv = gdf_csv.to_crs(epsg=out_epsg)
        gdf_crs.getFromEPSG(out_epsg)

    return gdf_csv, gdf_crs


def convert_csv_to_points(in_csvfilename, out_shapefilename=None, coord_columns=None,
                          coord_columns_epsg=4326, out_epsg=0):
    """ Converts a comma or tab delimited file to a shapefile.

    If EPSG is omitted it will default to 0 or unknown.

    Args:
        in_csvfilename (str):   Input CSV or Txt (tab) file.
        out_shapefilename (str):  Output Shapefile Name. If path is excluded it uses TEMPDIR
        coord_columns (list) : Column list representing coordinates. If None, it will be predicted.
        coord_columns_epsg (int):  EPSG number for specified coord_columns. Default is 4326 (WGS84)
        out_epsg (int):  The EPSG number representing the required output coordinates system.
                        - OR -
                        0  is no reprojection
                        -1 calculate utm zonal projection

    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.
        (pyprecag_crs): The coordinate system details of the dataframe

    Modified from:  https://gis.stackexchange.com/a/182234

    """
    start_time = time.time()

    if not os.path.exists(in_csvfilename):
        raise IOError("Invalid path: {}".format(in_csvfilename))

    for argCheck in [('coord_columns_epsg', coord_columns_epsg), ('out_epsg', out_epsg)]:
        if not isinstance(argCheck[1], (int, long)):
            raise TypeError('{} must be a integer.'.format(argCheck[0]))

    desc_csv = CsvDescribe(in_csvfilename)
    pdf_csv = desc_csv.open_pandas_dataframe()

    gdf_csv, gdf_crs = add_point_geometry_to_dataframe(pdf_csv, coord_columns,
                                                       coord_columns_epsg, out_epsg)

    if config.get_debug_mode() or out_shapefilename is not None:
        if out_shapefilename is None:
            out_shapefilename = '{}_{}.shp'.format(
                os.path.splitext(os.path.basename(in_csvfilename))[0],
                inspect.getframeinfo(inspect.currentframe())[2])

        save_geopandas_tofile(gdf_csv, out_shapefilename, overwrite=True)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=timedelta(seconds=time.time() - start_time)))

    return gdf_csv, gdf_crs


def convert_polygon_feature_to_raster(feature, pixel_size, value=1, nodata_val=-9999,
                                      snap_extent_to_pixel=True):
    """ Convert a feature into a raster (block grid)

    The "value" argument can be a numeric, or the name of a feature column. If the feature column
    is a string, then the index value will be used to assign a value to the raster.

    By default, snap_extent_to_pixel is true and is used to snap the bounding coordinates to
    a pixel edge.

    Args:
       feature (pandas.core.series.Series): a polygon feature to convert to a raster
       pixel_size (int, float]: the pixel size for the raster
       value (int, float, str): pixel value. A column from the feature can be used
       nodata_val (int, float): nodata value
       snap_extent_to_pixel (bool): round the extent coordinates to be divisible by the pixel size

    Returns:
       burned (numpy.ndarray) : The array representing the raster data
       Dict[str, int]] : The rasterio meta kwargs required for writing the array to disk.
    """
    if not isinstance(feature, (pd.Series, geopandas.GeoSeries)):
        raise TypeError("Input feature should be a geopandas series")

    if not isinstance(pixel_size, (int, long, float)):
        raise TypeError('Pixel size must be an integer or floating number.')

    if isinstance(value, str) and 'value' not in list(feature.index):
        raise TypeError('Value string {} is not a column name')
    elif not isinstance(value, (int, long, float)):
        raise TypeError('Value should be a column name, or a number')

    transform, width, height, bbox = create_raster_transform(feature.geometry.bounds,
                                                             pixel_size, snap_extent_to_pixel)

    gdf_feat = GeoDataFrame(GeoSeries(feature.geometry),
                            geometry=GeoSeries(feature.geometry).geometry)

    # create an affine transformation matrix to associate the array to the coordinates.
    if value in gdf_feat.columns:
        if is_string_dtype(gdf_feat[value]):
            shapes = ((geom, val) for geom, val in zip(gdf_feat.geometry, (gdf_feat.index + 1)))
        else:
            shapes = ((geom, val) for geom, val in zip(gdf_feat.geometry, gdf_feat[value]))
    else:
        shapes = ((geom, value) for geom in gdf_feat.geometry)

    burned = features.rasterize(shapes=shapes, out_shape=(height, width), fill=nodata_val,
                                transform=transform)
    meta = {}
    meta.update(height=height, width=width,
                transform=transform,
                nodata=nodata_val,
                dtype=rasterio.dtypes.get_minimum_dtype(burned))

    return burned, meta
