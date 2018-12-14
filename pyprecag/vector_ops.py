import collections
import datetime
import functools
import inspect
import logging
import os
import time
import warnings
from tempfile import NamedTemporaryFile

import fiona
import geopandas
import numpy as np
import pyproj
import shapely
from fiona import collection as fionacoll
from fiona.crs import to_string
from osgeo import osr
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from . import crs as pyprecag_crs
from . import TEMPDIR, config
from .describe import VectorDescribe, save_geopandas_tofile
from .errors import GeometryError

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
# LOGGER.setLevel(logging.DEBUG)
# DEBUG = config.get_debug_mode()  # LOGGER.isEnabledFor(logging.DEBUG))


def thin_point_by_distance(point_geodataframe, point_crs, thin_distance_metres=1.0, out_filename=None):
    """ Thin points by a distance in metres. All points less than the set distance will be removed.

    The point_crs must be a Projected Coordinate System

    Args:
        point_geodataframe (geopandas.geodataframe.GeoDataFrame): The input points geodataframe
        point_crs (pyprecag_crs.crs): The Projected Spatial Reference System of the point_geodataframe
        thin_distance_metres (float):   A floating number representing the minimum distance between points
        out_filename (str): (Optional) The path and filename to save the result to

    Returns:
        geopandas.geodataframe.GeoDataFrame: The thinned geodataframe

    TODO: Replace with scipy.spatial methods which doesn't rely on the file being correctly sorted.
          In depth testing required.
    """

    for argCheck in [('thinDist_m', thin_distance_metres)]:
        if not isinstance(argCheck[1], (int, long, float)):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not isinstance(point_geodataframe, geopandas.GeoDataFrame):
        raise TypeError("Input Points should be a geopandas data frame")

    if 'POINT' not in ','.join(list(point_geodataframe.geom_type.unique())).upper():
        raise GeometryError('Invalid geometry. input shapefile should be point or multipoint')

    if not isinstance(point_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if not point_crs.srs.IsProjected():
        raise TypeError("Input data is not in a projected coordinate system")

    if out_filename is not None:
        if not os.path.exists(os.path.dirname(out_filename)):
            raise IOError('Output directory {} does not exist'.format(os.path.dirname(out_filename)))

    if out_filename is None or config.get_debug_mode():
        # get a unique name, it will create, open the file and delete when done after saving the variable
        with NamedTemporaryFile(prefix='{}_'.format(inspect.getframeinfo(inspect.currentframe())[2]),
                                suffix='.shp', dir=TEMPDIR) as new_file:
            out_filename = new_file.name
        # out_filename = '{}.shp'.format(inspect.getframeinfo(inspect.currentframe())[2])

    # Create a copy
    point_geodataframe = point_geodataframe.copy()

    # create a column and fill with NAN
    point_geodataframe['filter'] = np.nan

    if thin_distance_metres == 0:
        LOGGER.warning('A thin distance of Zero (0) was used. Thinning of input data not undertaken')
        return point_geodataframe

    thinDistCSunits = thin_distance_metres
    if not point_crs.srs.IsProjected():

        # use lower left corner of bnd box to convert metres to dd94 distance.
        thinDistCSunits = pyprecag_crs.distance_metres_to_dd(point_geodataframe.total_bounds[0],
                                                             point_geodataframe.total_bounds[1], thin_distance_metres)

        LOGGER.warning('\Input data uses a geographics coordinate system.\n'
                       'Reprojecting parameter thin distance of {}m to {} for {}'.format(thin_distance_metres,
                                                                                         thinDistCSunits,
                                                                                         point_geodataframe.crs))

    else:
        LOGGER.debug('\nFiltering by distance {}m ({} for {}) ------'.format(thin_distance_metres,
                                                                             thinDistCSunits, point_geodataframe.crs))

    # create the function to use with pandas apply

    def thin_by_distance(curRow, filter_text):
        global prevPt
        if prevPt == '':
            prevPt = curRow
        else:
            dist = curRow.distance(prevPt)
            if dist <= thinDistCSunits:
                return filter_text
            else:
                prevPt = curRow

    global prevPt
    prevPt = ''
    subset = point_geodataframe[point_geodataframe['filter'].isnull()].copy()
    drop_cols = [ea for ea in subset.columns.tolist() if ea not in ['geometry', 'FID']]
    subset.drop(drop_cols, axis=1, inplace=True)

    # Update/Add x,y coordinates to match coordinate system
    for argCheck in ['pointX', 'pointY']:
        if argCheck in point_geodataframe.columns:
            warnings.warn('{} already exists. Values will be updated'.format(argCheck))

    subset['pointX'] = subset.geometry.apply(lambda p: p.x)
    subset['pointY'] = subset.geometry.apply(lambda p: p.y)

    iloop = 0
    for sortBy in ['pointXY', 'pointX', 'pointY']:
        filterTime = time.time()

        if sortBy != 'pointXY':
            subset.sort_values(by=sortBy, ascending=True, inplace=True)

        filter_string = '{} ({}m)'.format(sortBy, thinDistCSunits)
        subset['filter'] = subset['geometry'].apply(lambda x: thin_by_distance(x, filter_string))

        stepTotal = len(subset[subset['filter'] == filter_string])

        if stepTotal > 0:
            iloop += 1
            subset.loc[subset['filter'] == filter_string, 'filter_inc'] = iloop

            # update(join/merge) the master filter column with results from thinning.
            point_geodataframe.loc[point_geodataframe.index.isin(subset.index), 'filter'] = subset['filter']
            point_geodataframe.loc[point_geodataframe.index.isin(subset.index), 'filter_inc'] = subset['filter_inc']

        subset = subset[subset['filter'].isnull()].copy()
        if len(subset) == 1:
            raise TypeError(
                "There are no features left after {}. Check the coordinate systems and try again".format(sortBy))

        LOGGER.info(
            '{:<30} {:>10,}   {dur:<15} {}'.format('Filter by distance - {}'.format(filter_string.replace('point', '')),
                                                  len(subset), 'del {} pts'.format(stepTotal),
                                                  dur=datetime.timedelta(seconds=time.time() - filterTime)))

    # set sort back to original row order
    # point_geodataframe.sort_index(axis=1, ascending=True, inplace=True)

    filterTime = time.time()
    if config.get_debug_mode():  # save with filter column.
        save_geopandas_tofile(point_geodataframe, out_filename, overwrite=True)
    elif out_filename is not None:
        save_geopandas_tofile(point_geodataframe.drop('filter', axis=1), out_filename, overwrite=True)

    return point_geodataframe


def move_or_copy_vector_file(in_filename, out_filename, keepInput=True, overwrite=True):
    """ Move, Rename or Copy a vector file from one name/location to another.

       This uses the GDAL CopyDataSource and DeleteDataSource methods.
       To Move/Rename - keepInput=False.   It will copy the file to new location then delete original.
       To Copy - keepInput=True.

       Args:
           in_filename (str): The Input filename
           out_filename (str): The output filename and location
           keepInput (bool): if False will mimic a move/rename by copying the file then deleting the input file
           overwrite (bool): If output exists it will be overwritten.

       """
    start_time = time.time()

    if not os.path.exists(in_filename):
        raise IOError("Invalid path: {}".format(in_filename))

    if os.path.splitext(in_filename)[-1] != os.path.splitext(out_filename)[-1]:
        raise TypeError("Input and output should be the same file type, ie both shapefiles")

    if in_filename == out_filename:
        raise TypeError("Input and output should have different file names")

    if os.path.exists(out_filename) and not overwrite:
        raise IOError("Output Exists and overwrite is false: {}".format(in_filename))

    vectDesc = VectorDescribe(in_filename)
    geoDF = vectDesc.open_geo_dataframe()
    geoDF.to_file(out_filename, crs_wkt=vectDesc.crs.crs_wkt)

    LOGGER.debug('Successfully renamed from \n   {} \n    to    \n   {}'.format(in_filename, out_filename))
    LOGGER.info('{} complete !!  Duration H:M:SS - {dur}'.format(inspect.currentframe().f_code.co_name,
                                                              dur=datetime.timedelta(seconds=time.time() - start_time)))


def explode_multi_part_features(in_shapefilename, out_shapefilename):
    """ Convert Multipart Vector Features to Single Part Features.

    Args:
        in_shapefilename (str): The input vector file, normally a shapefile
        out_shapefilename (str): the resulting shapefile

    Returns:None

    """
    start_time = time.time()

    if not os.path.exists(in_shapefilename):
        raise IOError("Invalid path: {}".format(in_shapefilename))

    if os.path.splitext(in_shapefilename)[-1] != os.path.splitext(out_shapefilename)[-1]:
        raise TypeError("Input and output should be the same file type, ie both shapefiles")

    if os.path.basename(in_shapefilename) == os.path.basename(out_shapefilename):
        raise TypeError("Input and output should have different file names")

    vectDesc = VectorDescribe(in_shapefilename)

    if 'MULTI' not in vectDesc.geometry_type:
        LOGGER.exception('Input shapefiles does not contain multi-part geometry')
        return

    # open the original shapefile
    with fionacoll(in_shapefilename, 'r') as source:
        # create the new shapefile
        with fionacoll(out_shapefilename, 'w', driver=source.driver, crs=vectDesc.epsg,
                       schema=source.schema) as output:
            # iterate the features of the original shapefile
            for elem in source:
                # iterate the list of geometries in one element (split the MultiLineString)
                for line in shape(elem['geometry']):
                    # write the line to the new shapefile
                    output.write({'geometry': mapping(line), 'properties': elem['properties']})

    LOGGER.info('{} complete !!  Duration H:M:SS - {dur}'.format(inspect.currentframe().f_code.co_name,
                                                             dur=datetime.timedelta(seconds=time.time() - start_time)))


def calculate_area_length_in_metres(in_filename, dissolve_overlap=True):
    """Calculate the total Area and Length in metres for an input polygon file.

    If the input file is geographic, then the feature will be 'projected' into wgs84 utm system

    Currently there is no Summarize by Column or update existing features .. implemented (see TODO).

    Args:
        in_filename (str): Input polygon file.
        dissolve_overlap (bool): Return values after Dissolve ALL geometries.

    Returns:
        list[area,length]: The Total Area and Length
    """
    start_time = time.time()

    if not os.path.exists(in_filename):
        raise IOError("Invalid path: {}".format(in_filename))

    with fiona.open(in_filename, 'r') as source:
        inSRS = osr.SpatialReference()
        inSRS.ImportFromProj4(to_string(source.crs))

        if not inSRS.IsProjected():
            zone, utmSRS, wgs84SRS = pyprecag_crs.getUTMfromWGS84(*source.bounds[:2])

            # Create the Project Transformation object
            project = functools.partial(pyproj.transform,
                                        pyproj.Proj(**source.crs),
                                        pyproj.Proj(utmSRS.ExportToProj4()))

        # Store Results in a dictionary
        resultsDict = collections.defaultdict(float)

        # Store projected polygons
        polygonsGeom = []

        # loop through polygons collecting areas and lengths
        for eaFeat in source:
            # get the geometry of the feature
            eaGeom = shape(eaFeat['geometry'])

            # if this is geographic the it will be dd
            resultsDict['origArea'] += eaGeom.area
            resultsDict['origLength'] += eaGeom.length

            # If inSRS is Geographic the reproject the feature to a projected system to get area in metres
            if not inSRS.IsProjected():
                project = functools.partial(pyproj.transform,
                                            pyproj.Proj(**source.crs),
                                            pyproj.Proj(utmSRS.ExportToProj4()))

                # Get the projected geometry
                eaGeom = shapely.ops.transform(project, eaGeom)

            # store numbers in metres
            resultsDict['Area_m'] += eaGeom.area
            resultsDict['Length_m'] += eaGeom.length

            # If overlaps are to be dissolved we need a projected geometry collection
            if dissolve_overlap:
                polygonsGeom.append(eaGeom)

        if dissolve_overlap:
            # Run the dissolve.
            polygondis = shapely.ops.unary_union(polygonsGeom)
            resultsDict['Area_m_NoOverlap'] += polygondis.area
            resultsDict['Length_m_NoOverlap'] += polygondis.length

    # get the final values
    area = resultsDict['Area_m']
    length = resultsDict['Length_m']

    if dissolve_overlap:
        area = resultsDict['Area_m_NoOverlap']
        length = resultsDict['Length_m_NoOverlap']

    LOGGER.debug('CRS:              {}'.format(to_string(source.crs)))
    if not inSRS.IsProjected():
        LOGGER.debug('\t{:<25} Area (dd): {:>20.5f}   Length (dd): {:>20.5f}'.format('No Overlap-',
                                                                                     resultsDict['origArea'],
                                                                                     resultsDict['origLength']))
        LOGGER.debug('Projected To CRS: {}'.format(utmSRS.GetAttrValue("PROJCS", 0)))

    LOGGER.debug('\t{:<25} Area (m) : {:>20.5f}   Length (m) : {:>20.5f}'.format('Total with Overlap-',
                                                                                 resultsDict['Area_m'],
                                                                                 resultsDict['Length_m']))
    if dissolve_overlap:
        LOGGER.debug('\t{:<25} Area (m) : {:>20.5f}   Length (m) : {:>20.5f}'.format('Total with No Overlap-',
                                                                                     resultsDict['Area_m_NoOverlap'],
                                                                                     resultsDict['Length_m_NoOverlap']))

    LOGGER.debug('{} complete !!  Duration H:M:SS - {dur}'.format(inspect.currentframe().f_code.co_name,
                                                            dur=datetime.timedelta(seconds=time.time() - start_time)))
    return area, length
