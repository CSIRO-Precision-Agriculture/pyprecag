import collections
from datetime import timedelta
import functools
import inspect
import logging
import os
import time
import warnings
from tempfile import NamedTemporaryFile
import six

import fiona
from fiona import collection as fionacoll
from fiona.crs import to_string

import geopandas as gpd
import pyproj
from osgeo import osr

from ._compat import SHAPELY_GE_20
import shapely
from shapely import force_3d, remove_repeated_points
from shapely.geometry import mapping, shape

if SHAPELY_GE_20:
    from shapely import MultiPoint, LineString
else:
    from shapely.geometry import MultiPoint, LineString

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
        geopandas.GeoDataFrame: The thinned GeoDataFrame

    TODO: Replace with scipy.spatial methods which doesn't rely on the file being correctly sorted.
          In depth testing required.
    """

    for argCheck in [('thinDist_m', thin_distance_metres)]:
        if not isinstance(argCheck[1], six.integer_types + (float, )):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not isinstance(point_geodataframe, gpd.GeoDataFrame):
        raise TypeError("Input Points should be a geopandas data frame")

    if not any("POINT" in g.upper() for g in point_geodataframe.geom_type.unique()):
        raise GeometryError('Invalid geometry. input shapefile should be point or multipoint')

    if point_crs and not isinstance(point_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if (point_crs and not point_crs.srs.IsProjected()) or point_geodataframe.crs.is_geographic:
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
    # point_geodataframe['filter'] = np.nan

    if thin_distance_metres == 0:
        LOGGER.warning('A thin distance of Zero (0) was used. Thinning of input data not undertaken')
        return point_geodataframe

    thinDistCSunits = thin_distance_metres
    if point_geodataframe.crs.is_geographic:
        # use lower left corner of bnd box to convert metres to dd94 distance.
        thinDistCSunits = pyprecag_crs.distance_metres_to_dd(point_geodataframe.total_bounds[0],
                                                             point_geodataframe.total_bounds[1], thin_distance_metres)

        LOGGER.warning(f'Input data uses a geographics coordinate system. \nReprojecting parameter thin distance of '
                       f'{thin_distance_metres}m to {thinDistCSunits} for {point_geodataframe.crs}')

    else:
        LOGGER.debug(f'\nFiltering by distance {thin_distance_metres}m ({thinDistCSunits} for '
                     f'{point_geodataframe.crs}) ------')

    # Update/Add x,y coordinates to match coordinate system
    for argCheck in ['pointX', 'pointY']:
        if argCheck in point_geodataframe.columns:
            warnings.warn(f'{argCheck} already exists. Values will be updated')

    point_geodataframe['thinID'] = point_geodataframe.index

    subset = point_geodataframe.copy()
    subset.drop(columns=subset.columns.difference(['geometry', 'thinID']), axis=1, inplace=True)

    # convert to a 3d point using the thinID for z
    subset['geometry'] = force_3d(subset.geometry, subset['thinID'])
    iloop = 0
    subset['pointX'] = subset.geometry.x
    subset['pointY'] = subset.geometry.y
    subset.set_index('thinID', inplace=True)

    for sortBy in ['pointXY', 'pointX', 'pointY']:

        if sortBy != 'pointXY':
            subset.sort_values(by=sortBy, ascending=True, inplace=True)

        # play dot to dot to create a single very long line
        line = LineString(subset.geometry.tolist())

        # thin the line by the desired distance
        simp_line = remove_repeated_points(line, tolerance=thin_distance_metres)

        # turn vertex back into points
        multipt = MultiPoint(list(simp_line.coords))

        # extract the Z (or index) from the geometry
        thin_idx = [int(pt.z) for pt in multipt.geoms]

        filtered = subset.loc[~subset.index.isin(thin_idx)]

        LOGGER.info('remaining: {:.>10,} ... removed: {:.>10,} ... {}'.format(len(thin_idx), len(filtered),
                                                                              f'{sortBy} ({thin_distance_metres}m)'))
        if len(thin_idx) > 0:
            iloop += 1
            # update(join/merge) the master filter column with results from thinning.
            point_geodataframe.loc[point_geodataframe['thinID'].isin(filtered.index),
                        ['filter','filter_inc']] = f'{sortBy} ({thin_distance_metres}m)', iloop+1
        else:
            raise TypeError("There are no features left after {}. Check the coordinate systems and try again".format(sortBy))

        subset = subset.loc[subset.index.isin(thin_idx)].copy()


    # set sort back to original row order
    # point_geodataframe.sort_index(axis=1, ascending=True, inplace=True)

    if config.get_debug_mode():  # save with filter column.
        save_geopandas_tofile(point_geodataframe, out_filename, overwrite=True)
    elif out_filename is not None:
        save_geopandas_tofile(point_geodataframe.drop('filter', axis=1), out_filename, overwrite=True)
    point_geodataframe.drop(columns='thinID',inplace=True)
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
                                                              dur=str(timedelta(seconds=time.time() - start_time))))


def explode_multi_part_features(in_shapefilename, out_shapefilename):
    """ Convert Multipart Vector Features to Single Part Features.

    Args:
        in_shapefilename (str): The input vector file, normally a shapefile
        out_shapefilename (str): the resulting shapefile

    Returns:None

    """
    warnings.warn('explode_multi_part_features() is deprecated in favor of `GeoDataFrame.explode()`  '
                  'see https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.explode.html',
                  DeprecationWarning, stacklevel=2)

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
                                                             dur=str(timedelta(seconds=time.time() - start_time))))


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
    warnings.warn('calculate_area_length_in_metres() is deprecated in favor of '
                  '`gdf.to_crs(gdf.estimate_utm_crs().to_epsg()).geometry.area` or similar',
                  DeprecationWarning, stacklevel=2)

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
                                                            dur=(timedelta(seconds=time.time() - start_time))))
    return area, length
