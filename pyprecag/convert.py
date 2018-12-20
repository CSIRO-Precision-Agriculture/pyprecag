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
from . import general, describe, TEMPDIR, config
from describe import CsvDescribe, predictCoordinateColumnNames, VectorDescribe, save_geopandas_tofile
from .errors import GeometryError
from .raster_ops import RasterSnapExtent, create_raster_transform

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
DEBUG = config.get_config_key('debug_mode')  # LOGGER.isEnabledFor(logging.DEBUG)
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

def convertPolyToGrid(in_shapefilename,
                      out_rasterfilename,
                      pixel_size,
                      nodata_val=-9999,
                      snap_extent_to_pixel=True,
                      overwrite=True):
    """ Convert a polygon shapefile to a raster for a set pixel size.

    source : http://ceholden.github.io/open-geo-tutorial/python/chapter_4_vector.html

    Args:
        in_shapefilename (str): A polygon shapefile to rasterize
        out_rasterfilename (str): the output raster name
        pixel_size (float): the pixel size for the raster
        nodata_val (int):the value to use for nodata
        snap_extent_to_pixel (bool):round the extent coordinates to be divisible by the pixel size
        overwrite (bool): if true overwrite existing output file

    Returns:
        None:
    """

    start_time = time.time()
    # Add Function + arguments to Log
    if config.get_config_key('debug_mode'):
        [LOGGER.debug(ea) for ea in general.print_functions_string(inspect.currentframe())]

    if not os.path.exists(in_shapefilename):
        raise IOError("Invalid path: {}".format(in_shapefilename))

    # Open the data source and read in the extent
    shpDS = ogr.Open(in_shapefilename)
    shpLayer = shpDS.GetLayer()

    # Create the destination data source
    xMin, xMax, yMin, yMax = shpLayer.GetExtent()

    LOGGER.info('\tVector Extent: LL:{} {}'.format(xMin, yMin))
    LOGGER.info('\t               UR:{} {}'.format(xMax, yMax))

    # We may want to snap the output grids to a multiple of the grid size, allowing adjacent blocks to align nicely.
    if snap_extent_to_pixel:
        xMin, yMin, xMax, yMax = RasterSnapExtent(xMin, yMin, xMax, yMax, pixel_size)
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
        # raise Exception("Error rasterizing layer: %s" % err)

    else:
        LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                                   dur=datetime.timedelta(seconds=time.time() - start_time)))

    # close the data sources
    rBand = None
    rasterDS = None
    shpLayer = None
    shpDS = None

    return


def convertGridToVesper(in_rasterfilename, out_vesperfilename):
    """ Convert a Grid the format compatible with vesper.

       The vesper format contains the xy coordinates for the center of each pixel with one space before the x
       and three spaces between the x and y values.

       The coordinates are written as rows from the Upper Left corner of the raster.

       Args:
           in_rasterfilename (str):An input raster file
           out_vesperfilename (str): Save the output vesper file to this file.

    """
    start_time = time.time()
    # Add Function + arguments to Log
    if config.get_config_key('debug_mode'):
        [LOGGER.debug(ea) for ea in general.print_functions_string(inspect.currentframe())]

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

    # noinspection PyStringFormat
    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return


def convertPtToLines_fiona(in_filename, out_filename, thin_dist_m=1.0, aggregate_dist_m=25.0):
    """Convert a points vector file to lines, by connect-the-dot.

    This is used to create vehicle tracks from GPS data therefore the point order should be sorted by increasing
    time sequence to ensure the 'dot-to-dot' occurs correctly.

    For efficiency, points will be thinned. This will remove points less than the set thinDist. Resulting points will be
    connected to form lines.

    Line ends are detected where the distance between points is greater than the aggregate distance.

    Args:
        in_filename (str): Input vector points file normally a shapefile
        out_filename (str):  An output vector file of the same format as the in file.
        thin_dist_m (float): A floating number representing the minimum distance between points
        aggregate_dist_m (float): A floating number representing the maximum distance between point. This us used to
                        detect a line end.

    TODO: Replace point thinning alg with scipy.spatial method.

    """
    start_time = time.time()
    # Add Function + arguments to Log
    if config.get_config_key('debug_mode'):
        [LOGGER.debug(ea) for ea in general.print_functions_string(inspect.currentframe())]

    for argCheck in [('thin_dist_m', thin_dist_m), ('aggregate_dist_m', aggregate_dist_m)]:
        if not isinstance(argCheck[1], (int, long, float)):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not os.path.exists(in_filename):
        raise IOError("Invalid path: {}".format(in_filename))

    if os.path.splitext(in_filename)[-1] != os.path.splitext(out_filename)[-1]:
        raise TypeError("Input and output should be the same filetype, ie both shapefiles")

    if os.path.basename(in_filename) == os.path.basename(out_filename):
        raise TypeError("Input and output should have different filenames")

    if out_filename is None and config.get_config_key('debug_mode'):
        out_filename = '{}_{}.shp'.format(os.path.basename(in_filename)[:10],
                                          inspect.getframeinfo(inspect.currentframe())[2])

    # if out_shapefilename doesn't include a path then add tempdir.
    if out_filename is not None and not os.path.isabs(out_filename):
        out_filename = os.path.join(TEMPDIR, out_filename)

    vectDesc = VectorDescribe(in_filename)

    if 'POINT' not in vectDesc.geometry_type.upper():
        raise GeometryError('Invalid geometry. Input shapefile should be point or multipoint')

    thinDistCSunits = thin_dist_m
    aggregateDistCSunits = aggregate_dist_m

    # if input file is in geographics (assume WGS84) then reproject distances to coordinate system units
    if not vectDesc.crs.srs.IsProjected():
        thinDistCSunits = pyprecag_crs.distance_metres_to_dd(vectDesc.extent[0], vectDesc.extent[1], thin_dist_m)
        aggregateDistCSunits = pyprecag_crs.distance_metres_to_dd(vectDesc.extent[0], vectDesc.extent[1], aggregate_dist_m)

    schema = {'geometry': 'LineString', 'properties': {'id': 'int'}}
    ptCount = 0
    with fionacoll(in_filename, "r") as pointInput:
        # Using crs=pointInput.crs incorrectly applies the coordinate system as sutm or gda2020 so
        # use the describe value instead
        with fionacoll(out_filename, "w", driver=pointInput.driver, schema=schema,
                       crs=from_epsg(vectDesc.crs.epsg)) as lineOutput:
            line = ()
            lines = []
            lineAdded = False  # has current line been added to multiline
            lineCount = 0
            for curPoint in pointInput:
                if len(line) == 0:
                    line += ((curPoint['geometry']['coordinates']),)  # Add the first Point
                else:

                    ptDist = shape(curPoint['geometry']).distance(shape(prevPoint['geometry']))

                    if ptDist > thinDistCSunits:
                        if ptDist <= aggregateDistCSunits:
                            line += ((curPoint['geometry']['coordinates']),)
                            lineAdded = False
                        else:  # start a new line
                            if len(line) > 1:
                                lineCount += 1
                                # lines.append(line)
                                ptCount += len(line)
                                lineAdded = True
                                lineOutput.write({'properties': {'id': lineCount},
                                                  'geometry': mapping(LineString(line))})

                            line = ((curPoint['geometry']['coordinates']),)  # start a fresh collection

                prevPoint = curPoint

            # add the last line
            if not lineAdded and len(line) > 1:
                lineCount += 1
                lineOutput.write({'properties': {'id': lineCount},
                                  'geometry': mapping(LineString(line))})

        if lineCount == 0:
            raise GeometryError(
                "Empty Geometry: The output file contains no features. Check the coordinate system and try again")
            # LOGGER.exception("Empty Geometry - the output file contains no features. "
            #                  "Check the Coordinate System and try again",True)

        LOGGER.debug('Created {} lines using {} points from a total of {} points\n'
                     'Using a thinning distance of {} m  and an aggregate distance of {} m '.format(
            lineCount, ptCount, vectDesc.feature_count, thin_dist_m, aggregate_dist_m))

        # noinspection PyStringFormat
        LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                                   dur=datetime.timedelta(seconds=time.time() - start_time)))

    return out_filename


def addPointGeometryToDataframe(in_dataframe, coord_columns=None, coord_columns_EPSG=0, out_EPSG=0):
    """ Add Point geometry to a pandas or geopandas data frame through sets of coordinate columns.
       If required, it will also reproject to the chosen EPGS coordinate system

    Args:
        in_dataframe (DataFrame): The input dataframe. This could be a pandas or geopandas dataframe
        coord_columns (list) : A list of the columns representing coordinates. If None, it will be predicted.
        coord_columns_EPSG (int):  The EPSG number representing the coordinates inside the csv file
        out_EPSG (int): The EPSG number representing the required output coordinates system.
                        Use  0 for no projection
    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.

    Modified from:  Make a Shapefile from a geopandas dataframe - https://gis.stackexchange.com/a/182234

    """
    start_time = time.time()
    # Add Function + arguments to Log
    if config.get_config_key('debug_mode'):
        [LOGGER.debug(ea) for ea in general.print_functions_string(inspect.currentframe())]

    if not isinstance(in_dataframe, (geopandas.GeoDataFrame, pd.DataFrame)):
        raise TypeError("Input points should be a geopandas dataframe")

    for argCheck in [('coord_columns_EPSG', coord_columns_EPSG), ('out_EPSG', out_EPSG)]:
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
    gdfCSV = geopandas.GeoDataFrame(in_dataframe, geometry='geometry', crs=from_epsg(coord_columns_EPSG))

    # drop the original geometry columns to avoid confusion
    gdfCSV.drop([x_column, y_column], axis=1, inplace=True)

    gdfCRS = pyprecag_crs.crs()
    gdfCRS.getFromEPSG(coord_columns_EPSG)

    if out_EPSG == -1:
        xmin, ymin, _, _ = gdfCSV.total_bounds
        out_EPSG =  pyprecag_crs.getProjectedCRSForXY(xmin, ymin, coord_columns_EPSG).epsg_number

    if out_EPSG > 0:
        gdfCSV = gdfCSV.to_crs(epsg=out_EPSG)
        gdfCRS.getFromEPSG(out_EPSG)

    return gdfCSV, gdfCRS


def convertCsvToPoints(in_csvfilename, out_shapefilename=None, coord_columns=None, coord_columns_EPSG=4326, out_EPSG=0):
    r""" Converts a comma or tab delimited file to a shapefile.

    If EPSG is omitted it will default to 0 or unknown.

    Args:
        in_csvfilename (str):   Input CSV or Txt (tab) file.
        out_shapefilename (str):  Output Shapefile Name. If the path is excluded it will save to TEMPDIR
        coord_columns (list) : A list of the columns representing coordinates. If None, it will be predicted.
        coord_columns_EPSG (int):  The EPSG number for the specified coord_columns. Default is 4326 (WGS84)
        out_EPSG (int):  The EPSG number representing the required output coordinates system.
                        - OR -
                        0  is no reprojection
                        -1 calculate utm zonal projection

    Returns:
        (GeoDataFrame): A geodataframe of Points with ALL attributes.
        (pyprecag_crs): The coordinate system details of the dataframe

    Modified from:  Make a Shapefile from a geopandas dataframe - https://gis.stackexchange.com/a/182234
    Examples:
        >>> gpf, crs = convertCsvToPoints(r"../test/data/area2_yield_file_ISO-8859-1.csv",coord_columns_EPSG=4326)
        >>> type(gpf)
        <class 'geopandas.geodataframe.GeoDataFrame'>
        >>> len(gpf)
        10000
        >>> gpf.total_bounds
        array([ 142.35507623,  -35.65127353,  142.36843437,  -35.64060432])

    """
    start_time = time.time()
    # Add Function + arguments to Log
    if config.get_config_key('debug_mode'):
        [LOGGER.debug(ea) for ea in general.print_functions_string(inspect.currentframe())]

    if not os.path.exists(in_csvfilename):
        raise IOError("Invalid path: {}".format(in_csvfilename))

    for argCheck in [('coord_columns_EPSG', coord_columns_EPSG), ('out_EPSG', out_EPSG)]:
        if not isinstance(argCheck[1], (int, long)):
            raise TypeError('{} must be a integer.'.format(argCheck[0]))

    descCSV = CsvDescribe(in_csvfilename)
    pdfCSV = descCSV.open_pandas_dataframe()

    gdfCSV, gdfCRS = addPointGeometryToDataframe(pdfCSV, coord_columns, coord_columns_EPSG, out_EPSG)

    if config.get_config_key('debug_mode') or out_shapefilename is not None:
        if out_shapefilename is None:
            out_shapefilename = '{}_{}.shp'.format(os.path.splitext(os.path.basename(in_csvfilename))[0],
                                                   inspect.getframeinfo(inspect.currentframe())[2])
        save_geopandas_tofile(gdfCSV, out_shapefilename, overwrite=True)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return gdfCSV, gdfCRS


def convertPolygonFeatureToRaster(feature,  pixel_size, value=1, nodata_val=-9999,
                                  snap_extent_to_pixel=True,
                                  overwrite=True):

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
       overwrite (bool):if true overwrite existing output file

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

    transform, width, height, bbox = create_raster_transform(feature.geometry.bounds,pixel_size)

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
