# -*- coding: utf-8 -*-
"""
Format conversion routines.
"""
import datetime
import inspect
import logging
import os
import time

import geopandas
import pandas as pd
import rasterio
from fiona import collection as fionacoll
from fiona.crs import from_epsg
from geopandas import GeoDataFrame, GeoSeries
from osgeo import ogr, gdal
from pandas.core.dtypes.common import is_string_dtype
from rasterio import features
from shapely.geometry import Point, mapping, shape, LineString

from . import crs as pyprecag_crs
from . import describe, TEMPDIR, config
from .config import DEBUG
from describe import CsvDescribe, predictCoordinateColumnNames, VectorDescribe, save_geopandas_tofile
from .errors import GeometryError
from .raster_ops import raster_snap_extent, create_raster_transform

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
# DEBUG = config.get_config_key('debug_mode')  # LOGGER.isEnabledFor(logging.DEBUG)
# LOGGER.setLevel("DEBUG")


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

    if pixelsize_metres % 1000.0 == 0:   # or km
        result = '{}km'.format(int(pixelsize_metres / 1000))
    elif pixelsize_metres % 1 == 0:     # is in m
        result = '{}m'.format(int(pixelsize_metres))
    elif (pixelsize_metres * 100.0) % 1 == 0:     # is in cm
        result = '{}cm'.format(int(pixelsize_metres * 100))
    else:  # or mm
        result = '{}mm'.format(int(pixelsize_metres * 1000))
    return result


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

    shpDS = ogr.Open(in_shapefilename)
    shpLayer = shpDS.GetLayer()

    xMin, xMax, yMin, yMax = shpLayer.GetExtent()

    LOGGER.info('\tVector Extent: LL:{} {}'.format(xMin, yMin))
    LOGGER.info('\t               UR:{} {}'.format(xMax, yMax))

    # We may want to snap the output grids to a multiple of the grid size, allowing adjacent blocks to align nicely.
    if snap_extent_to_pixel:
        xMin, yMin, xMax, yMax = raster_snap_extent(xMin, yMin, xMax, yMax, pixel_size)
        LOGGER.info('\tRaster Extent: Snap is {}'.format(snap_extent_to_pixel))
        LOGGER.info('\t           LL:{} {}'.format(xMin, yMin))
        LOGGER.info('\t           UR:{} {}'.format(xMax, yMax))

    # This '+1' is not a standard approach, but gives compatible rows and cols with arcpy
    xCols = int((xMax - xMin) / pixel_size) + 1
    yRows = int((yMax - yMin) / pixel_size) + 1
    LOGGER.info('\t           Cols: {}'.format(xCols))
    LOGGER.info('\t           Rows: {}'.format(yRows))

    # for pixel data types see http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4
    tifDrv = gdal.GetDriverByName('GTiff')
    if tifDrv is None:
        # If the above Fails then here is an alternate version.
        tifDrv = describe.Find_GDALDriverByName('GTiff')

    rasterDS = tifDrv.Create(out_rasterfilename, xCols, yRows, 1, gdal.GDT_Int16)

    # Define spatial reference
    rasterDS.SetProjection(shpLayer.GetSpatialRef().ExportToWkt())
    rasYMax = yMin + (yRows * pixel_size)
    rasterDS.SetGeoTransform((xMin, pixel_size, 0, rasYMax, 0, -pixel_size))

    rBand = rasterDS.GetRasterBand(1)
    rBand.SetNoDataValue(nodata_val)
    rBand.Fill(nodata_val)

    # Rasterize the shapefile layer to our new dataset
    # options = ['ALL_TOUCHED=TRUE','ATTRIBUTE=id']
    err = gdal.RasterizeLayer(rasterDS,  # output to our new dataset
                              [1],  # output to our new datasets first band
                              shpLayer,  # rasterize this layer
                              None, None,  # don't worry about transformations since we're in same projection
                              burn_values=[1],  # burn value 1
                              options=['ALL_TOUCHED=FALSE'])

    if err != 0:
        LOGGER.exception("Error rasterizing layer: {}".format(err), True)

    else:
        LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                                   dur=datetime.timedelta(seconds=time.time() - start_time)))

    # close the data sources
    rBand = None
    rasterDS = None
    shpLayer = None
    shpDS = None

    return


def convert_grid_to_vesper(in_rasterfilename, out_vesperfilename):
    """ Convert a Grid the format compatible with vesper.

       The vesper format contains the xy coordinates for the center of each pixel with one space before the x
       and three spaces between the x and y values.

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
                        xy = oRasterDS.xy(iRowY, iColX)  # takes rows,cols and returns an xy coord for cell center
                        vesFile.write(' {}   {}\n'.format(xy[0], xy[1]))

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return

def add_point_geometry_to_dataframe(in_dataframe, coord_columns=None, coord_columns_epsg=0, out_epsg=0):
    """ Add Point geometry to a pandas or geopandas data frame through sets of coordinate columns.
       If required, it will also reproject to the chosen EPSG coordinate system

    Args:
        in_dataframe (DataFrame): The input dataframe. This could be a pandas or geopandas dataframe
        coord_columns (list) : A list of the columns representing coordinates. If None, it will be predicted.
        coord_columns_epsg (int):  The EPSG number representing the coordinates inside the csv file
        out_epsg (int): The EPSG number representing the required output coordinates system.
                        Use  0 for no projection
    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.

    Modified from:  Make a Shapefile from a geopandas dataframe - https://gis.stackexchange.com/a/182234

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
    in_dataframe['geometry'] = in_dataframe.apply(lambda x: Point((float(x[x_column]), float(x[y_column]))), axis=1)

    # Now, convert the pandas DataFrame into a GeoDataFrame. The geopandas constructor expects a geometry column which
    # can consist of shapely geometry objects, so the column we created is just fine:
    gdfCSV = geopandas.GeoDataFrame(in_dataframe, geometry='geometry', crs=from_epsg(coord_columns_epsg))

    # drop the original geometry columns to avoid confusion
    gdfCSV.drop([x_column, y_column], axis=1, inplace=True)

    gdfCRS = pyprecag_crs.crs()
    gdfCRS.getFromEPSG(coord_columns_epsg)

    if out_epsg == -1:
        xmin, ymin, _, _ = gdfCSV.total_bounds
        out_epsg =  pyprecag_crs.getProjectedCRSForXY(xmin, ymin, coord_columns_epsg).epsg_number

    if out_epsg > 0:
        gdfCSV = gdfCSV.to_crs(epsg=out_epsg)
        gdfCRS.getFromEPSG(out_epsg)

    return gdfCSV, gdfCRS


def convert_csv_to_points(in_csvfilename, out_shapefilename=None, coord_columns=None, coord_columns_epsg=4326, out_epsg=0):
    """ Converts a comma or tab delimited file to a shapefile.

    If EPSG is omitted it will default to 0 or unknown.

    Args:
        in_csvfilename (str):   Input CSV or Txt (tab) file.
        out_shapefilename (str):  Output Shapefile Name. If the path is excluded it will save to TEMPDIR
        coord_columns (list) : A list of the columns representing coordinates. If None, it will be predicted.
        coord_columns_epsg (int):  The EPSG number for the specified coord_columns. Default is 4326 (WGS84)
        out_epsg (int):  The EPSG number representing the required output coordinates system.
                        - OR -
                        0  is no reprojection
                        -1 calculate utm zonal projection

    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.
        (pyprecag_crs): The coordinate system details of the dataframe

    Modified from:  Make a Shapefile from a geopandas dataframe - https://gis.stackexchange.com/a/182234

    """
    start_time = time.time()

    if not os.path.exists(in_csvfilename):
        raise IOError("Invalid path: {}".format(in_csvfilename))

    for argCheck in [('coord_columns_epsg', coord_columns_epsg), ('out_epsg', out_epsg)]:
        if not isinstance(argCheck[1], (int, long)):
            raise TypeError('{} must be a integer.'.format(argCheck[0]))

    descCSV = CsvDescribe(in_csvfilename)
    pdfCSV = descCSV.open_pandas_dataframe()

    gdfCSV, gdfCRS = add_point_geometry_to_dataframe(pdfCSV, coord_columns, coord_columns_epsg, out_epsg)

    if DEBUG or out_shapefilename is not None:
        if out_shapefilename is None:
            out_shapefilename = '{}_{}.shp'.format(os.path.splitext(os.path.basename(in_csvfilename))[0],
                                                   inspect.getframeinfo(inspect.currentframe())[2])
        save_geopandas_tofile(gdfCSV, out_shapefilename, overwrite=True)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return gdfCSV, gdfCRS


def convert_polygon_feature_to_raster(feature, pixel_size, value=1, nodata_val=-9999, snap_extent_to_pixel=True):

    """ Convert a feature into a raster (block grid)

    The "value" argument can be a numeric, or the name of a feature column. If the feature column is a string, then the
    index value will be used to assign a value to the raster.

    By default, snap_extent_to_pixel is true and is used to snap the bounding coordinates to a pixel edge.

    Args:
       feature (pandas.core.series.Series): a polygon feature to convert to a raster
       pixel_size (int, float]: the pixel size for the raster
       value (int, float, str): the value to assign to the pixel. A column from the feature can be used
       nodata_val (int, float): the value to use for nodata
       snap_extent_to_pixel (bool): round the extent coordinates to be divisible by the pixel size

    Returns:
       burned (numpy.ndarray) : The array representing the raster data
       Dict[str, int]] : The rasterio meta kwargs required for writing the array to disk.
    """
    if not isinstance(feature, (pd.Series,geopandas.GeoSeries)):
        raise TypeError("Input feature should be a geopandas series")

    if not isinstance(pixel_size, (int, long, float)):
        raise TypeError('Pixel size must be an integer or floating number.')

    if isinstance(value,str) and 'value' not in list(feature.index):
        raise TypeError('Value string {} is not a column name')
    elif not isinstance(value,(int,long,float)):
        raise TypeError('Value should be a column name, or a number')

    transform, width, height, bbox = create_raster_transform(feature.geometry.bounds, pixel_size, snap_extent_to_pixel)

    gdfFeat = GeoDataFrame(GeoSeries(feature.geometry), geometry=GeoSeries(feature.geometry).geometry)

    # create an affine transformation matrix to associate the array to the coordinates.
    if value in gdfFeat.columns:
        if is_string_dtype(gdfFeat[value]):
            shapes = ((geom, val) for geom, val in zip(gdfFeat.geometry, (gdfFeat.index + 1)))
        else:
            shapes = ((geom, val) for geom, val in zip(gdfFeat.geometry, gdfFeat[value]))
    else:
        shapes = ((geom, value) for geom in gdfFeat.geometry)

    burned = features.rasterize(shapes=shapes, out_shape=(height,width), fill=nodata_val, transform=transform)
    meta = {}
    meta.update(height=height, width=width,
                transform=transform,
                nodata=nodata_val,
                dtype=rasterio.dtypes.get_minimum_dtype(burned))

    return burned, meta
