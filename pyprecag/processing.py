from datetime import timedelta
import glob
import inspect
import re
from collections import defaultdict
import warnings
import logging
import os
from pathlib import Path

import six
from six.moves import zip as izip

import random
import time

from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame, GeoSeries
from osgeo import gdal

import rasterio
from rasterio import crs as riocrs
from rasterio import features
from rasterio.io import MemoryFile
from rasterio.fill import fillnodata
from rasterio.mask import mask as rio_mask
from rasterio.warp import reproject, Resampling
from rasterio.windows import get_data_window, intersection, from_bounds

from scipy.cluster.vq import *
from scipy import stats

from ._compat import SHAPELY_GE_20

if SHAPELY_GE_20:
    from shapely import Point, LineString, force_2d, MultiLineString
    from shapely.geometry import mapping
else:
    from shapely.geometry import Point, LineString, mapping
    from .convert import drop_z as force_2d

from shapely.ops import linemerge, substring

from . import TEMPDIR, config, crs as pyprecag_crs, number_types
from .crs import from_epsg
from .table_ops import calculate_strip_stats
from .convert import (convert_grid_to_vesper, numeric_pixelsize_to_string, convert_polygon_feature_to_raster,
                      deg_to_8_compass_pts, line_bearing, point_to_point_bearing, text_rotation)

from .raster_ops import (focal_statistics, save_in_memory_raster_to_file, reproject_image,
                         calculate_image_indices, stack_and_clip_rasters)

from .describe import save_geopandas_tofile, VectorDescribe, predictCoordinateColumnNames
from .errors import GeometryError
from .vector_ops import thin_point_by_distance

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def block_grid(in_shapefilename, pixel_size, out_rasterfilename, out_vesperfilename,
               out_epsg=0, groupby=None, nodata_val=-9999, snap=True, overwrite=False):
    """Convert a polygon boundary to a 0,1 raster and generate a VESPER compatible list of
       coordinates for kriging.

       The input polygon shapefile will be reproject to the out_epsg coordinate system.

       Will process in a loop based on the groupby field name. All polygons with the same
       attribute will be dissolved together and a single block grid will be created. The
       groupby attribute will be added to the output filename.

        Args:
            in_shapefilename (str): Input polygon shapefile
            pixel_size (float):  The required output pixel size
            out_rasterfilename (str): Filename of the raster Tiff that will be created
            out_vesperfilename (str):  The output vesper file
            out_epsg (int): The epsg number for the output raster coordinate system.
            groupby (str): Column name to group polygon features with.
            nodata_val (int): an integer to use as nodata
            snap (bool): Snap Extent to a factor of the Pixel size
            overwrite (bool): if true overwrite existing file

    """

    if not isinstance(pixel_size, number_types):
        raise TypeError('Pixel size must be an integer or floating number.')

    for ea_arg in [('nodata_val', nodata_val), ('out_epsg', out_epsg)]:
        if not isinstance(ea_arg[1], six.integer_types):
            raise TypeError('{} must be a integer - Got {}'.format(*ea_arg))

    if out_epsg <= 0:
        raise ValueError('EPSG Code ({}) must be a positive integer above 0'.format(out_epsg))

    if riocrs.CRS.from_epsg(out_epsg).is_geographic:
        raise ValueError('EPSG Code ({}) must be a projected coordinate system'.format(out_epsg))

    desc_poly_shp = VectorDescribe(in_shapefilename)

    gdf_poly = desc_poly_shp.open_geo_dataframe()

    if 'POLY' not in desc_poly_shp.geometry_type.upper():
        raise GeometryError('Invalid Geometry. Input shapefile should be polygon or multipolygon')

    if groupby is None or groupby.strip() == '':
        groupby = None

    if groupby is not None and groupby not in gdf_poly.columns:
        raise ValueError('Groupby column {} does not exist'.format(groupby))

    # -------------------------------------------------------------------------------------------
    # reproject shapefile
    if out_epsg != 0 and desc_poly_shp.crs.epsg_number != out_epsg:
        gdf_poly.to_crs(epsg=out_epsg, inplace=True)

    del desc_poly_shp

    # If a column name is specified, dissolve by this and create multi-polygons
    # otherwise dissolve without a column name at treat all features as one
    if groupby is not None:
        # change null/nones to blank string
        gdf_poly[groupby].fillna('', inplace=True)

        gdf_poly = gdf_poly.dissolve(groupby, as_index=False).copy()

    else:
        gdf_poly = GeoDataFrame(geometry=[gdf_poly.unary_union])

    # Loop through polygons features ----------------------------------------------------------------
    output_files = []
    for index, feat in gdf_poly.iterrows():
        loop_time = time.time()
        step_time = time.time()
        status = '{} of {}'.format(index + 1, len(gdf_poly))

        if groupby is not None:
            feat_name = re.sub('[^0-9a-zA-Z]+', '-', str(feat[groupby])).strip('-')
            if feat_name == '':
                feat_name = 'No-Name'

            r_file, r_ext = os.path.splitext(out_rasterfilename)
            out_file = r_file + "_" + feat_name + r_ext
        else:
            out_file = out_rasterfilename
            feat_name = 'All Polygons'

        # From shapes create a blockgrid mask. -----------------------------------------------------
        blockgrid, new_blockmeta = convert_polygon_feature_to_raster(feat, pixel_size, snap_extent_to_pixel=snap)

        with rasterio.open(out_file, 'w', driver='GTiff', count=1, tfw='YES',
                           crs=from_epsg(out_epsg), **new_blockmeta) as dest:
            dest.write(blockgrid, 1)
        output_files += [out_file]

        del blockgrid, new_blockmeta

        LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format('Feature to block_grid', status, feat_name,
                                                          dur=str(timedelta(seconds=time.time() - step_time))))

        r_file, r_ext = os.path.splitext(out_file)
        convert_grid_to_vesper(in_rasterfilename=out_file,
                               out_vesperfilename=out_file.replace(r_ext, '_v.txt'))

    return output_files


def create_polygon_from_point_trail(points_geodataframe, points_crs, out_filename, thin_dist_m=1.0,
                                    aggregate_dist_m=25, buffer_dist_m=10, shrink_dist_m=3):
    """Create a polygon from a Point Trail created from a file containing GPS coordinates.

    The point order should be sorted by increasing time sequence to ensure the 'dot-to-dot' occurs
    correctly.

    The workflow is as follows:
          Points -> Lines -> Buffer Out (expand) -> Buffer In (shrink)

    For efficiency, points will be thinned. This will remove points less than the set thinDist.
    Resulting points will be connected to form lines.

    Line ends are detected where the distance between points is greater than the aggregate distance.
     This will occur at a turning point or interruption (ie creek) and gps point collection is
    stopped. Typically this distance is slightly greater than the row/swath width.

    Buffering is used to convert to polygon and remove the gap between row_count/swaths. The buffer
    distance is usually half the row/swath width.

    Buffering with a negative value is used to remove excess area on the outside of the polygon.
    This shrink distance is usually around 7m less than the buffer distance.

    Args:
        points_geodataframe (geopandas.geodataframe.GeoDataFrame): Input points vector geodataframe
        points_crs (pyprecag_crs.crs): The Projected Spatial Reference System of the
                    point_geodataframe
        out_filename (str): Output polygon file. This should be the same format as the input.
        thin_dist_m (float): The minimum distance in metres between points to be used to thin the
                    points dataset.
        aggregate_dist_m (int): A floating number representing the maximum distance between point.
                    This is used to detect a line end. Typically this is slightly larger than
                    the row/swath width.
        buffer_dist_m (int): The Buffer distance in metres. Typically half the swath or row width.
        shrink_dist_m (int): The shrink distance in metres. Typically about 7 less than the buffer
                    distance.
    """
    warnings.warn('points_crs as parameter and return values are deprecated in favor of `geopandas.crs` and '
                  'will be removed in a future version', FutureWarning, stacklevel=2)

    for argCheck in [('thin_dist_m', thin_dist_m), ('aggregate_dist_m', aggregate_dist_m),
                     ('buffer_dist_m', buffer_dist_m), ('shrink_dist_m', shrink_dist_m)]:
        if not isinstance(argCheck[1], number_types):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not isinstance(points_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : inputGeoDataFrame')

    if not any("POINT" in g.upper() for g in points_geodataframe.geom_type.unique()):
        raise TypeError('Invalid input data : A points geopandas dataframe is required')

    if  points_crs and not isinstance(points_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if out_filename is None or out_filename == '':
        raise TypeError('Please specify an output filename')

    if not os.path.exists(os.path.dirname(out_filename)):
        raise IOError('Output directory {} does not exist'.format(os.path.dirname(out_filename)))

    # if points_crs:
    #     points_geodataframe.crs = points_crs.epsg

    start_time = time.time()
    points_geodataframe = points_geodataframe.copy()
    if 'FID' not in points_geodataframe.columns:
        points_geodataframe['FID'] = points_geodataframe.index

    # don't need any attribution so drop it all except fid and geometry
    dropcols = [ea for ea in points_geodataframe.columns.tolist() if ea not in ['geometry', 'FID']]
    points_geodataframe.drop(dropcols, axis=1, inplace=True)

    gdf_thin = thin_point_by_distance(points_geodataframe, None, thin_dist_m)
    gdf_thin = gdf_thin[gdf_thin['filter'].isnull()].copy()

    del points_geodataframe

    dropcols = [ea for ea in gdf_thin.columns.tolist() if ea not in ['geometry', 'FID']]
    gdf_thin.drop(dropcols, axis=1, inplace=True)

    step_time = time.time()
    # calc distance between two points
    gdf_thin['dist_shift'] = gdf_thin.distance(gdf_thin.shift(1))
    # Calc first value to be 0
    gdf_thin.loc[gdf_thin['dist_shift'].isnull(), 'dist_shift'] = 0

    # apply a new row/line id when distance between points is greater than 25m
    gdf_thin["lineID"] = (gdf_thin["dist_shift"] >= aggregate_dist_m).cumsum()
    LOGGER.info('{:<30} {:<15} {dur}'.format('Assign lineID', '',
                                             dur=str(timedelta(seconds=time.time() - step_time))))

    temp_file_list = []
    if config.get_debug_mode():
        temp_file_list = [os.path.join(TEMPDIR, os.path.basename(
            out_filename.replace('.shp', '_0thnpts.shp')))]

        save_geopandas_tofile(gdf_thin, temp_file_list[-1], overwrite=True)

    step_time = time.time()
    # A line must have at least two points so remove lines with only one point.
    # count all points per lineID
    ptsperline = gdf_thin.groupby('lineID').count()

    # select all LineID's with more than one point
    gdf_thin = gdf_thin.loc[gdf_thin['lineID'].isin(ptsperline[ptsperline['geometry'] > 1].index)]

    # convert to lines
    # source https://gis.stackexchange.com/a/202274
    try:
        # Aggregate these points with the GroupBy - for shapely 1.3.2+ use
        df_line = gdf_thin.groupby(['lineID'])['geometry'].apply(lambda x: LineString(x.tolist()))
    except ValueError:
        df_line = gdf_thin.groupby(['lineID'])['geometry'].apply(
            lambda x: LineString([(p.x, p.y) for p in x]))

    gdf_final = GeoDataFrame(df_line, geometry='geometry',crs=gdf_thin.crs)
    del gdf_thin, ptsperline, df_line

    LOGGER.info('{:<30} {:<15} {dur}'.format('Convert to lines', '',
                                             dur=str(timedelta(seconds=time.time() - step_time))))

    if config.get_debug_mode():
        temp_file_list.append(os.path.join(TEMPDIR, os.path.basename(
            out_filename.replace('.shp', '_1line.shp'))))

        save_geopandas_tofile(gdf_final, temp_file_list[-1], overwrite=True)

    # buffer and dissolve overlap of lines
    step_time = time.time()
    gdf_final = GeoDataFrame(geometry=[gdf_final.buffer(buffer_dist_m).unary_union],crs=gdf_final.crs)
    # gdf_final.crs = ptsgdf_crs
    gdf_final['FID'] = gdf_final.index
    LOGGER.info('{:<30} {:<15} {dur}'.format('Buffer by {}'.format(buffer_dist_m), '',
                                             dur=str(timedelta(seconds=time.time() - step_time))))

    if config.get_debug_mode():
        temp_file_list.append(os.path.join(
            TEMPDIR, os.path.basename(out_filename.replace('.shp', '_2buf.shp'))))

        save_geopandas_tofile(gdf_final, temp_file_list[-1], overwrite=True)

    if shrink_dist_m != 0:
        step_time = time.time()
        gdf_final = GeoDataFrame(geometry=gdf_final.buffer(-abs(shrink_dist_m)),crs=gdf_final.crs)

        LOGGER.info('{:<30} {:<15} {dur}'.format('Shrink by {}'.format(shrink_dist_m), '',
                                                 dur=str(timedelta(seconds=time.time() - step_time))))

        if config.get_debug_mode():
            temp_file_list.append(os.path.join(
                TEMPDIR, os.path.basename(out_filename.replace('.shp', '_3bufshrink.shp'))))
            save_geopandas_tofile(gdf_final, temp_file_list[-1], overwrite=True)

    step_time = time.time()

    gdf_final = GeoDataFrame(geometry=gdf_final.geometry.explode(index_parts=False),crs=gdf_final.crs)

    gdf_final['FID'] = gdf_final.index
    gdf_final['Area'] = gdf_final.area
    gdf_final['Perimeter'] = gdf_final.length

    LOGGER.info('{:<30} {:<15} {dur}'.format('Explode polygons', '',
                                             dur=str(timedelta(seconds=time.time() - step_time))))

    save_geopandas_tofile(gdf_final, out_filename, overwrite=True)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=str(timedelta(seconds=time.time() - start_time))))

    thin_ratio = (4 * 3.14 * gdf_final['Area'].sum() /
                  (gdf_final['Perimeter'].sum() * gdf_final['Perimeter'].sum()))

    LOGGER.debug('{:<25} {:.5f}'.format('Thin Ratio : ', thin_ratio))

    if thin_ratio < 0.01:
        LOGGER.warning('For an improved result, increase the buffer width and shrink '
                       'distance and try again.')
        return 'For an improved result, increase the buffer width and shrink distance and ' \
               'try again.'
    else:
        return gdf_final


def clean_trim_points(points_geodataframe, points_crs, process_column, output_csvfile, poly_geodataframe=None,
                      boundary_polyfile=None, out_keep_shapefile=None, out_removed_shapefile=None,
                      remove_zeros=True, stdevs=3, iterative=True, thin_dist_m=1.0):
    """ Clean and/or Trim a points dataframe.

        Preparation includes:
            - Clip data to polygon.
            - Move data to a Projected Coordinate system.
            - Remove values where data_column are less than or equal to zero
            - Calculate Normalised value for data_column (number of StDev).
            - Iteratively Trim outliers based on Normalised data_column
            - Remove points closer than a set distance (trim_dist_m)

    Args:
        points_geodataframe (geopandas.geodataframe.GeoDataFrame): The input points geodataframe
        points_crs (pyprecag_crs.crs): The Spatial Reference System of the point_geodataframe
        process_column (str):  The column  to normalise, trim and clean.
        output_csvfile (str): The Trimmed & Cleaned output CSV file
        out_keep_shapefile (str): Optionally save the Trimmed & Cleaned to a shapefile
        out_removed_shapefile (str): Optionally save a shapefile containing features removed while
                    cleaning/filtering. A column called filter will be added showing the reason a
                    point was removed.
        boundary_polyfile (str): Optionally a polygon used to Clip the points.
        remove_zeros (bool): Optionally remove values where data_column are <= to zero
                    prior to removing outliers.
        stdevs (int): The number of standard deviations used to trim outliers
        iterative (bool): Optionally Iteratively Trim outliers based on Normalised Data column
        thin_dist_m (float): A distance in metres representing the minimum allowed distance
                    between points. Points less than this distance will be removed.

    Returns:
        geopandas.geodataframe.GeoDataFrame: Representing the cleaned/trimmed data file.
                   New columns included in the output are:
                      Eastings, Northings  - The point coordinates in a projected coordinate system
                      EN_Epsg  - The epsg number of the projected coordinate system
                      Nrm_    - normalise data_column calculation.
        pyprecag_crs.crs: The pyprecag CRS object of the points dataframe.

    """
    warnings.warn('points_crs as parameter and return values are deprecated in favor of `geopandas.crs` and '
                  'will be removed in a future version', FutureWarning, stacklevel=2)

    if boundary_polyfile is not None:
        warnings.warn('boundary_polyfile will be removed in future release. Please use poly_geodataframe instead',
                      FutureWarning, stacklevel=2)

    if not isinstance(points_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : points_geodataframe')

    if not any("POINT" in g.upper() for g in points_geodataframe.geom_type.unique()):
        raise TypeError('Invalid input data : a points geopandas dataframe is required')

    if points_crs and not isinstance(points_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if output_csvfile is None:
        raise TypeError('Invalid input data : an output csv path and file is required')

    for argCheck in [output_csvfile, out_keep_shapefile, out_removed_shapefile]:
        if argCheck is None:
            continue

        if not os.path.exists(os.path.dirname(argCheck)):
            raise IOError('Output folder does not exist: {}'.format(os.path.dirname(argCheck)))

    for argCheck in [('thinDist_m', thin_dist_m), ('stdev', stdevs)]:
        if not isinstance(argCheck[1], number_types):
            raise TypeError('{} must be a integer or floating number.'.format(argCheck[0]))

    if process_column not in points_geodataframe.columns:
        raise TypeError('Column {} does not exist in input geodataframe'.format(process_column))

    for argCheck in [('remove_zeros', remove_zeros), ('iterative', iterative)]:
        if not isinstance(argCheck[1], bool):
            raise TypeError('{} should be a boolean.'.format(argCheck[0]))

    if points_crs:
        points_geodataframe.crs = points_crs.epsg

    norm_column = 'nrm_' + process_column
    LOGGER.info('Normalized Column is {}'.format(norm_column))

    if norm_column in points_geodataframe.columns:
        LOGGER.warning('Column {} already exists and will be overwritten'.format(norm_column))

    if boundary_polyfile is not None:
        warnings.warn('boundary_polyfile will be removed in future release. Please use poly_geodataframe instead',
                      FutureWarning, stacklevel=2)

        if not os.path.exists(boundary_polyfile):
            raise IOError("Invalid path: {}".format(boundary_polyfile))
        poly_geodataframe = GeoDataFrame.from_file(boundary_polyfile)

    if poly_geodataframe is not None and not isinstance(poly_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : poly_geodataframe')

    if not any("POLY" in g.upper() for g in poly_geodataframe.geom_type.unique()):
        raise GeometryError('Invalid geometry. Input poly_geodataframe or boundary_polyfile should be polygon'
                            ' or multipolygon')

    start_time = time.time()

    # set a UniqueID Field which ISNT the FID for use through out the processing
    id_col = 'PT_UID'  # IE clean/trim fid
    points_geodataframe[id_col] = points_geodataframe.index
    gdf_points = points_geodataframe.copy()

    # To speed things up, drop all un-required columns
    gdf_points.drop(columns=gdf_points.columns.difference(['geometry', id_col, norm_column, process_column]),
                    axis=1, inplace=True)

    step_time = time.time()

    # set defaults
    gdf_points['filter'] = np.nan
    gdf_points['filter_inc'] = 0

    def add_filter_message(filter_string):
        # use ..... to space out values as QGIS doesn't honour multiple spaces or tabs in the log panel.
        if filter_string in gdf_points['filter'].unique():
            total_count = len(gdf_points)
            del_count = gdf_points.fillna('').value_counts(['filter_inc', 'filter'],sort=False).iloc[-1]
            left_count = gdf_points.fillna('').value_counts(['filter_inc', 'filter'],sort=False).iloc[0]
            # left_count = len(gdf_points[gdf_points['filter'].isnull()])
            # del_count = len(gdf_points[gdf_points['filter'] == filter_string])

            LOGGER.info(f'remaining: {left_count:.>10,} ... removed: {del_count:.>10,} ... {filter_string}')

    # Remove rows where data col is empty/null
    gdf_points.loc[gdf_points[process_column].isnull(),
                        ['filter', 'filter_inc']] = ['null/missing data', len(gdf_points['filter'].unique())]

    add_filter_message('null/missing data')

    # remove missing or empty geometry
    gdf_points.loc[(gdf_points.is_empty | gdf_points.geometry.isna()), ['filter', 'filter_inc']] = [
        'empty/missing geom', len(gdf_points['filter'].unique())]
    add_filter_message('empty/missing geom')

    # Remove duplicated geometries
    # https://github.com/geopandas/geopandas/issues/521#issuecomment-382806444
    gdf_points.loc[(gdf_points['filter'].isnull()) & (gdf_points["geometry"].to_wkb().duplicated(keep='first')),
                        ['filter', 'filter_inc']] = ['Duplicate XY', len(gdf_points['filter'].unique())]

    add_filter_message('Duplicate XY')

    if poly_geodataframe is not None:
        if poly_geodataframe.crs.to_epsg() != gdf_points.crs.to_epsg():
            # Preference is to the projected coordinate system then only project the smaller
            # dataset (usually poly) By now points should be in the out projected coordinate system
            poly_geodataframe.to_crs(gdf_points.crs, inplace=True)

            LOGGER.info(f'{"Reproject clip polygon": <30} {"to": >10}   {gdf_points.crs.to_epsg():<15}')
            # dur=str(timedelta(seconds=time.time() - step_time))))

        step_time = time.time()

        # Clip to boundary then apply to filter column
        gdf_points.loc[(gdf_points['filter'].isnull()) & (~gdf_points.geometry.within(poly_geodataframe.unary_union)),
                            ['filter', 'filter_inc']] = ['clip', len(gdf_points['filter'].unique())]
        add_filter_message('clip')

        if gdf_points['filter'].isnull().sum() == 0:
            raise GeometryError('Clipping removed all features. Check coordinate systems and/or '
                                'clip polygon layer and try again')

        del poly_geodataframe

        step_time = time.time()

    if remove_zeros:
        gdf_points.loc[(gdf_points['filter'].isnull()) & (gdf_points[process_column] <= 0),
                            ['filter', 'filter_inc']] = ['<= zero', len(gdf_points['filter'].unique())]

        add_filter_message('<= zero')

        if gdf_points['filter'].isnull().sum() == 0:
            raise GeometryError(f"Zero filter removed all points in column {process_column}")

    if stdevs > 0:
        i = 0
        # Get a fresh copy of subset.
        subset = gdf_points.loc[gdf_points['filter'].isnull()].copy()
        while len(subset) > 0:
            i += 1

            subset = gdf_points.loc[gdf_points['filter'].isnull()].copy()
            filter_str = f'{stdevs} std iter {int(i)}'

            yld_mean = subset[process_column].mean()
            yld_std = subset[process_column].std()

            subset[norm_column] = subset[process_column].apply(lambda x: (x - yld_mean) / yld_std)
            subset = subset.loc[subset[norm_column].abs() >= stdevs]

            gdf_points.loc[subset.index, ['filter', 'filter_inc']] = [filter_str, len(gdf_points['filter'].unique())]

            add_filter_message(filter_str)

            if not iterative:
                break

    gdf_thin = thin_point_by_distance(gdf_points[gdf_points['filter'].isnull()], None, thin_dist_m)
    # update the filter incremental number
    gdf_thin['filter_inc'] = gdf_thin['filter_inc'] + len(gdf_points['filter'].dropna().unique())

    # update(join/merge) gdfPoints['filter'] column with results from thinning.
    gdf_points.update(gdf_thin[gdf_thin['filter'].notnull()])
    step_time = time.time()

    del gdf_thin

    # Add the incremental number to the filter column ie '01 clip'
    gdf_points['filter'] = gdf_points[['filter_inc', 'filter']].dropna().apply(
        lambda x: '{:02n} {}'.format(x[0], x[1]), axis=1)

    # recalculate normalised col.
    yld_mean = gdf_points[gdf_points['filter'].isnull()][process_column].mean()
    yld_std = gdf_points[gdf_points['filter'].isnull()][process_column].std()
    gdf_points.loc[gdf_points['filter'].isnull(), norm_column] = (gdf_points[process_column] - yld_mean) / yld_std

    # prepare some summary results for filtered features.
    # Filter is the reason a point is removed,and filter_inc keeps them in the order they were
    # removed. ie std it before thining by distance.
    results_table = gdf_points.copy()
    results_table.loc[gdf_points['filter'].isnull(), ['filter_inc', 'filter']] = \
        [len(gdf_points['filter'].unique()), 'Pts remaining']

    results_table = results_table[['filter_inc', 'filter']].value_counts(sort=False).to_frame('count')
    results_table.reset_index(drop=False, inplace=True)

    results_table.loc[results_table['filter'].isnull(),
        ['filter_inc', 'filter']] = [len(results_table), 'Pts remaining']

    results_table['%'] = ((results_table['count'] / results_table['count'].sum()) * 100).round(3)
    results_table.sort_values('filter_inc', inplace=True)

    results_table.loc['Total'] = results_table.sum(numeric_only=True)  # this will convert int's to floats !!
    results_table.loc[results_table['filter'].isnull(), ['filter_inc', 'filter']] = [len(results_table), 'Total']

    # Clean up filtered results by removing all columns except those new ones which
    # have to be copied back to original

    dropcols = [ea for ea in gdf_points.columns.tolist()
                if ea not in [ norm_column, 'filter', id_col]]

    gdf_points.drop(dropcols, axis=1, inplace=True)

    # Clean up the original input dataframe and remove existing geometry and coord columns
    alt_coord_columns = config.get_config_key('geoCSV')['xCoordinate_ColumnName']
    alt_coord_columns += config.get_config_key('geoCSV')['yCoordinate_ColumnName']

    # Find and Drop coord columns if already exist.
    coord_columns = [fld for fld in points_geodataframe.columns if fld.upper() in alt_coord_columns]
    # coord_columns = coord_columns + ['geometry']

    if len(coord_columns) > 0:
        points_geodataframe.drop(coord_columns, axis=1, inplace=True)

    # Use geopandas merge instead of concat to maintain coordinate system info etc.
    gdf_final = points_geodataframe.merge(gdf_points, on=id_col, how='left')
    # Add x,y coordinates to match coordinate system
    gdf_final['Easting'] = gdf_final.geometry.x
    gdf_final['Northing'] = gdf_final.geometry.y
    gdf_final['EN_EPSG'] = gdf_final.crs.to_epsg()

    # move newer columns to end (geometry, filter, etc)
    gdf_final.drop([id_col], inplace=True, axis=1)

    # get the appropriate file encoding and save the file
    # TODO: CHECK THIS I dont think setting a default in the config will work.
    try:
        file_encoding = config.read_config()['geoCSV']['file_encoding']
    except KeyError:
        file_encoding = 'utf-8'
    # file_encoding = get_dataframe_encoding(gdf_final)

    gdf_final[gdf_final['filter'].isnull()].drop(['geometry', 'filter'], axis=1) \
        .to_csv(output_csvfile, index=False, encoding=file_encoding)

    LOGGER.info(f'{"Save to CSV":<30} {len(gdf_final[gdf_final["filter"].isnull()]): >10,} '
                f'  {os.path.basename(output_csvfile):<15}'
                f' {str(timedelta(seconds=time.time() - step_time))}')

    step_time = time.time()

    # gdfPoints have all points with a filter string assigned
    if out_keep_shapefile is not None and out_keep_shapefile != '':
        save_geopandas_tofile(gdf_final[gdf_final['filter'].isnull()].drop(['filter'], axis=1),
                              out_keep_shapefile, overwrite=True)

        LOGGER.info('{:<30} {:>10,}   {:<15} {dur}'.format(
            'Save to Shapefile', len(gdf_final[gdf_final['filter'].isnull()]),
            os.path.basename(out_keep_shapefile), dur=str(timedelta(seconds=time.time() - step_time))))

        step_time = time.time()

    if out_removed_shapefile is not None and out_removed_shapefile != '':
        if gdf_final[gdf_final['filter'].notnull()].empty:
            LOGGER.info('No features removed. Shapefile containing removed points not created.')
        else:
            save_geopandas_tofile(gdf_final[gdf_final['filter'].notnull()]
                                  .drop([norm_column], axis=1),
                                  out_removed_shapefile, overwrite=True)

            LOGGER.info('{:<30} {:>10,}   {:<15} {dur}'
                        .format('Save removed Shapefile', len(gdf_final[gdf_final['filter'].notnull()]),
                                os.path.basename(out_removed_shapefile),
                                dur=str(timedelta(seconds=time.time() - step_time))))

    LOGGER.info('\nResults:---------------------------------------\n{}\n'.format(
        results_table.drop('filter_inc', axis=1).to_string(index=False, justify='center',
                                                           formatters={'filter_inc': "{:,.0f}".format,
                                                                       'count'     : "{:,.0f}".format,
                                                                       '%'         : "{:.3f}%".format})))

    LOGGER.info(f'{process_column} mean....{yld_mean:.5f} ')
    LOGGER.info(f'{process_column} std....{yld_std: .5f} ')
    LOGGER.info(f'{process_column} CV....{ 100 * yld_std / yld_mean: .5f} ')

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=str(timedelta(seconds=time.time() - start_time))) )

    warnings.warn('return value points_crs is deprecated. Please use  `pyproj, geopandas solutions`',
                  FutureWarning, stacklevel=2)

    if points_crs :
        points_crs = pyprecag_crs.crs()
        points_crs.getFromEPSG(gdf_final.crs.to_epsg())

    return gdf_final[gdf_final['filter'].isnull()], points_crs


def random_pixel_selection(raster, raster_crs, num_points, out_shapefile=None):
    """Select randomly distributed valid data pixels from a raster file and convert to points
    representing the center of the pixel. There is an option to save to shapefile if required.

    Note:
        When opening a raster using RasterIO the coordinate system is imported from the Proj4 string
         and converted to a crs_wkt. This means that the crs_wkt recorded against the opened dataset
          does not always equal the crs_wkt of the original raster file. To remedy this, use
          pyprecag_crs.getCRSfromRasterFile to create a crs object

    Args:
        raster (rasterio.io.DatasetReader):Opened raster file via rasterio.open(os.path.normpath())
        num_points (int): The number of random sample points to select.
        raster_crs (pyprecag_crs.crs): The Spatial Reference System for the raster file
        out_shapefile (str):Optional.. the path and name of a shapefile used to save the points.

    Returns:
        geopandas.geodataframe.GeoDataFrame: A dataframe containing the select pixels as points
        pyprecag_crs.crs: The pyprecag CRS object of the points dataframe.

    """

    if not isinstance(raster, rasterio.DatasetReader):
        raise TypeError("Input should be a rasterio.DatasetReader created "
                        "using rasterio.open(os.path.normpath())")

    if not isinstance(num_points, six.integer_types):
        raise TypeError('Size must be an Integer.')

    if not isinstance(raster_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if out_shapefile is not None:
        if not os.path.exists(os.path.dirname(out_shapefile)):
            raise IOError('Output directory {} does not exist'.format(
                os.path.dirname(out_shapefile)))

    elif config.get_debug_mode():
        out_shapefile = '{}_randompts.shp'.format(inspect.getframeinfo(inspect.currentframe())[2])

    band1 = raster.read(1, masked=True)

    if num_points > band1.count():
        raise ValueError('Size parameter is greater than the number of data pixels')

    # Sources:  enumerate using a masked array https://stackoverflow.com/a/8621145
    # Create a geodataframe while iterating.... (Cell 2)
    #    http://nbviewer.jupyter.org/github/geopandas/geopandas/blob/master/examples/overlays.ipynb

    # Create a flattened mask array
    mask = ~band1.mask.ravel()

    # create an array with the row/col address and current pixel value for data pixels only
    # np.ndenumerate(image_read)   returns ((row,col),val) which can be unpacked
    # using ((r,c),v) and converted to a list using (int(r),int(c),v)

    band1_index = [(int(r), int(c), v) for ((r, c), v), m in izip(np.ndenumerate(band1), mask) if m]

    # get a fixes random sample
    rand_samp = random.sample(band1_index, num_points)

    # add coordinates for each sample point
    rand_samp_xy = [pt + raster.xy(pt[0], pt[1]) for pt in rand_samp]

    # convert to geopandas dataframe. Drop the data value from the output.
    random_pts_gdf = GeoDataFrame([{'geometry': Point(pt[3], pt[4]), 'PtID': i, 'X': pt[3],
                                    'Y'       : pt[4]}
                                   for i, pt in enumerate(rand_samp_xy)], crs=raster_crs.epsg)

    if out_shapefile is not None:
        save_geopandas_tofile(random_pts_gdf, out_shapefile, overwrite=True)

    return random_pts_gdf, raster_crs


def extract_pixel_statistics_for_points(points_geodataframe, points_crs, rasterfiles,
                                        output_csvfile, function_list=[np.nanmean], size_list=[3]):
    """Extract statistics from a list of rasters at set locations.

    All raster files in the list should be of the same pixel size. While multi-bands raster files
    are supported as an input, statistics will only be calculated and extracted for the first band.

    Statistics are calculated on pixel values with a square neighbourhood and saved to a CSV file.

    Pixels assigned with Nodata are converted to np.nan and excluded from the calculations.

    All original columns are included in the output csv file in addition to columns containing the
    values for each raster -> size -> statistic combination.  The column names can either be
    specified or derived via the out_colname argument in raster_ops.focal_statistics.

    A size of 1 can be used to extract exact pixel values and whereby no statistics are calculated.

    Args:
        points_geodataframe (geopandas.geodataframe.GeoDataFrame): The input points geodataframe of
                   locations to extract statistics.
        points_crs (pyprecag_crs.crs): The Spatial Reference System of the point_geodataframe
        rasterfiles (List[str]): the list of paths & file names for the input rasters
        output_csvfile (str): the path and filename of the output CSV.
        function_list (List[function]): A list of statistical functions to apply to the raster.
                   These can include numpy functions like np.nanmean or custom ones like pixel_count
        size_list (List[int]): The list of neighbourhood sizes used to apply statistical filtering.

    Returns:
        geopandas.geodataframe.GeoDataFrame: dataframe of the points and calculated statistics
        pyprecag_crs.crs: The pyprecag CRS object of the points dataframe.
    """
    warnings.warn('points_crs as parameter and return values are deprecated in favor of `geopandas.crs` and '
                  'will be removed in a future version', FutureWarning, stacklevel=2)

    if not isinstance(points_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : inputGeodataFrame')

    if not any("POINT" in g.upper() for g in points_geodataframe.geom_type.unique()):
        raise TypeError('Invalid input data : a points geopandas dataframe is required')

    if points_crs and not isinstance(points_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if output_csvfile is None or output_csvfile == '':
        raise TypeError('Please specify an output CSV filename')

    if not os.path.isabs(output_csvfile):
        # if a path isn't provided then create in temp. Used for daisy chaining functions.
        output_csvfile = os.path.join(TEMPDIR, output_csvfile)

    if not os.path.exists(os.path.dirname(output_csvfile)):
        raise IOError('Output directory {} does not exist'.format(os.path.dirname(output_csvfile)))

    if not isinstance(function_list, list):
        raise TypeError('Invalid Type: functions should be a list of function_list '
                        'statistics to process')

    for ea_function in function_list:
        # check to see if the functions in the list are functions ie are callable not random
        # string or numbers....
        if not callable(ea_function):
            raise TypeError('Invalid Type: should be a numpy or custom function like np.nanstd etc')

    if not isinstance(size_list, list):
        raise TypeError('Invalid Type: size_list should be a list of sizes to process')

    for ea_size in size_list:
        if not isinstance(ea_size, int) or ea_size % 2 == 0:
            raise TypeError("Size {} should be an odd integer. Only square "
                            "filters are supported.".format(ea_size))

    if not isinstance(rasterfiles, list):
        raise TypeError('Invalid Type: rasterfiles should be a list')

    if isinstance(rasterfiles[0], tuple):
        not_exists = [my_file for my_file, name in rasterfiles if not os.path.exists(my_file)]
    else:
        not_exists = [my_file for my_file in rasterfiles if not os.path.exists(my_file)]

    if len(not_exists) > 0:
        raise IOError('rasterfiles: {} raster file(s) do not exist\n\t({})'.format(
            len(not_exists), '\n\t'.join(not_exists)))

    # Check that all raster files have the same pixel size_list.
    pix = None
    res_error = []
    for ea_raster in rasterfiles:
        if isinstance(ea_raster, tuple):
            ea_raster, _ = ea_raster

        with rasterio.open(ea_raster) as src:
            if pix is None:
                pix = src.res
            elif pix != src.res:
                res_error.append(ea_raster)

    if len(res_error) > 0:
        raise TypeError('rasterfiles: Inconsistent pixel sizes'
                        '\n\t{}'.format('\n\t'.join(res_error)))

    '''TODO: consider clipping raster to extent of the points file which will be quicker when
            the raster is larger than the points.'''

    start_time = time.time()

    # drop null geometry
    points_geodataframe.dropna(subset=['geometry'], axis=0, inplace=True)
    original_crs = points_geodataframe.crs
    # overwrite the gdf proj4 string with the epsg mapping equivalent
    # if points_crs:
    #      points_geodataframe.crs = points_crs

    cur_pixel = False
    if 1 in size_list:
        cur_pixel = True
        size_list.remove(1)

    for ea_raster in rasterfiles:

        if isinstance(ea_raster, tuple):
            ea_raster, prefix = ea_raster
        else:
            prefix = ''

        # Need to get the wktproj of the raster from gdal NOT rasterio.
        # RasterIO works from the proj4 string NOT the wkt string so aussie zones details gets lost.
        rast_crs = pyprecag_crs.getCRSfromRasterFile(ea_raster)

        if rast_crs.epsg != points_geodataframe.crs.to_epsg():
            # we have to reproject the vector points to match the raster
            print('WARNING: Projecting points from {} to {}'.format(points_geodataframe.crs,
                                                                    rast_crs.epsg))
            points_geodataframe.to_crs(epsg=rast_crs.epsg_number, inplace=True)

        focal_stats = []
        with rasterio.open(ea_raster) as src:
            meta = src.meta.copy()
            if cur_pixel:
                focal_stats += [focal_statistics(src, size=1, out_colname=prefix)]

            for ea_size in size_list:
                for ea_function in function_list:
                    focal_stats += [focal_statistics(src, size=ea_size, function=ea_function, out_colname=prefix)]

        # Update the metadata for the desired focal_stats outputs.
        meta.update({'driver': 'GTiff', 'nodata': -9999, 'crs': rast_crs.crs_wkt,
                     'dtype' : rasterio.float32, 'count': len(focal_stats)})

        # create a raster in memory
        with MemoryFile() as memfile:
            # open the in-memory raster for writing
            with memfile.open(**meta) as dest:
                # loop through the outputs and write bands and tags
                for i, (arr, nm) in enumerate(focal_stats, start=1):
                    dest.write(arr, i)
                    dest.update_tags(i, name=nm)

            ''' https://gis.stackexchange.com/a/190428
             Get a generator object of  XY, value pairs to extract which can be used later by
              rasterio.sample() to extract values from ALL bands .... '''
            pt_shapes = ((geom.x, geom.y) for geom in points_geodataframe.geometry)

            # open the file for reading
            with memfile.open() as dest:
                # using rasterio sample, extract values from each band at each point
                raster_vals = [list(val) for val in dest.sample(pt_shapes)]

                # extract the statistic type from the band tag name element
                col_names = [dest.tags(i_band)['name'] for i_band in range(1, dest.count + 1)]

            del pt_shapes

        # Convert raster_vals numpy array to dataframe using tags as column names
        df_raster_vals = pd.DataFrame(raster_vals, columns=col_names)

        # replace values for points outside the raster extent with np.nan values
        df_raster_vals.replace(meta['nodata'], np.nan, inplace=True)

        del raster_vals

        # Add the point geometry back in by joining by row index
        # If Rounding's required do this.
        # points_geodataframe = pd.concat([points_geodataframe, df_raster_vals.round(6)], axis=1)
        points_geodataframe = pd.concat([points_geodataframe, df_raster_vals], axis=1)
        del df_raster_vals

    # Reproject points back to original coordinate system if required.
    if original_crs != points_geodataframe.crs:
        points_geodataframe.to_crs(original_crs, inplace=True)

    # Make sure the output CSV contains coordinates
    if None in predictCoordinateColumnNames(points_geodataframe.columns.tolist()):
        # probably from a shapefile so add coords. 
        if points_geodataframe.crs.is_geographic:
            points_geodataframe['Longitude'] = points_geodataframe.geometry.apply(lambda p: p.x)
            points_geodataframe['Latitude'] = points_geodataframe.geometry.apply(lambda p: p.y)

        else:
            points_geodataframe['Easting'] = points_geodataframe.geometry.apply(lambda p: p.x)
            points_geodataframe['Northing'] = points_geodataframe.geometry.apply(lambda p: p.y)

    # and for good measure add the associated epsg number
    points_geodataframe['EPSG'] = points_geodataframe.crs.to_epsg()

    if output_csvfile is not None:
        step_time = time.time()
        # Save to CSV only points to KEEP using appropriate file encoding
        points_geodataframe.drop(['geometry'], axis=1).to_csv(output_csvfile, index=False)

        LOGGER.info('{:<30}\t{:>10}   {dur:<15} {}'.format(
            'Saved CSV', '', os.path.basename(output_csvfile),
            dur=str(timedelta(seconds=time.time() - step_time))
        ))

    LOGGER.info('{:<30}\t{:>10} {dur:>15}\t{}'.format(
        inspect.currentframe().f_code.co_name, '', '',
        dur=str(timedelta(seconds=time.time() - start_time))
    ))


    if points_crs :
        points_crs = pyprecag_crs.crs()
        points_crs.getFromEPSG(points_geodataframe.crs.to_epsg())

    return points_geodataframe, points_crs


def multi_block_bands_processing(image_file, pixel_size, out_folder, band_nums=[], image_epsg=0,
                                 image_nodata=0, polygon_shapefile=None, groupby=None):
    """Derive multiple resampled image bands matching the specified pixel size and block grid extent
     for each shapefile polygon.

    Use this tool create individual band images for each polygon within a shapefile. A group-by
    column may be used to dissolve multiple polygons belonging to an individual block. The fitting
    of rasters to a base Block (grid) ensures for easier, more accurate multi-layered analysis
    required by in Precision Agriculture.

    The processing steps to achieve this are:
        - Dissolve polygons optionally using the group-by column

        Loop through each polygon feature and......
          - Create Block Grid
          - Clip image to polygon
          - Resample and fit to block grid for a pixel size and using Average resampling technique
          - Identify holes and fill if necessary
          - smooth image using a 5x5 pixel moving average (focal_statistics)

    If a polygon shapefile is not specified, polygons will be created from the images' mask

    The polygon shapefile will be re-projected to match the image file

    A block grid will be created for each feature for the nominated pixel size and used as the base
    for analysis

    "image_epsg" and "image_nodata" can be used to set the coordinate system and image nodata values
     when they are not present within the image file.

    The output filename will consist of the selected band number or the value of a band's rasterio
    custom name tag

    Args:
        image_file (str): An input image
        pixel_size (int, float): The desired output pixel size in metres.
        out_folder (str): The output folder for the created images.
        band_nums  (List[int]): a list of band numbers to process. If empty all bands will be used
        image_epsg (int): the epsg number for the image. Only used if not provided by the image.
        image_nodata (int): the nodata value for the image. Only used if not provided by the image.
        polygon_shapefile (str): a polygon shapefile used to cut up an image.
        groupby (str):  the column/field to use to group multiple features.

    Returns:
        List[str]: a list of created files.
    """

    if isinstance(polygon_shapefile, str) and polygon_shapefile.strip() == '':
        polygon_shapefile = None
    if isinstance(groupby, str) and groupby.strip() == '':
        groupby = None

    for ea_arg in [('image_file', image_file), ('out_folder', out_folder),
                   ('polygon_shapefile', polygon_shapefile)]:
        if ea_arg[0] != 'polygon_shapefile':
            if ea_arg[1] is None or ea_arg[1].strip() == '':
                raise ValueError('{} is required '.format(ea_arg[0]))

        if ea_arg[1] is not None and not os.path.exists(ea_arg[1]):
            raise IOError('{} does not exist-Got {}'.format(ea_arg[0], os.path.dirname(ea_arg[1])))

    for ea_arg in [('pixel_size', pixel_size), ('image_nodata', image_nodata)]:
        if not isinstance(ea_arg[1], number_types):
            raise TypeError('{} must be a integer or floating number - Got {}'.format(*ea_arg))

    if not isinstance(image_epsg, six.integer_types):
        raise TypeError('image_epsg must be a integer - Got {}'.format(image_epsg))

    if not isinstance(band_nums, list):
        raise TypeError('band_nums must be a list of bands to resample')

    if len(band_nums) == 0:
        band_nums = list(rasterio.open(image_file).indexes)
    else:
        image_bands = list(rasterio.open(image_file).indexes)
        band_check = [val for val in band_nums if val not in image_bands]
        if len(band_check) > 0:
            raise ValueError('{} are invalid bands.'.format(', '.join(band_check)))

    if polygon_shapefile is not None:
        vect_desc = VectorDescribe(polygon_shapefile)

        if 'POLYGON' not in vect_desc.geometry_type.upper():
            raise GeometryError('Invalid geometry. Input shapefile should be poly or multipolygon')

        if groupby is not None and groupby not in vect_desc.get_column_names():
            raise ValueError('Input groupby column does not exist in the shapefile')

        del vect_desc

    pixel_size_str = numeric_pixelsize_to_string(pixel_size)
    filename, ext = os.path.splitext(os.path.basename(image_file))
    temp_file_list = []
    # ----------------------------------------------------------------------------------------------
    # Run coordinate system checks on the input image
    with rasterio.open(os.path.normpath(image_file)) as src:
        if src.crs is None:
            if image_epsg is None or image_epsg == 0:
                raise ValueError('Input coordinate system required - image_file does not contain'
                                 'a coordinate system, and in_epsg is 0')
        else:
            image_epsg = pyprecag_crs.getCRSfromRasterFile(image_file).epsg_number

    # ----------------------------------------------------------------------------------------------
    # run checks on input polygon, and if required reproject it to the raster coordinate system
    gdf_poly = None
    if groupby is None or groupby.strip() == '':
        groupby = None

    if polygon_shapefile is None:
        # create a polygon from the image. hopefully by now the nodata val is correct.
        with rasterio.open(image_file) as src:
            # source: https://gis.stackexchange.com/a/187883
            # Read dataset's valid data mask as a ndarray and extract from all bands (dataset_mask).
            band = src.dataset_mask()
            mask = band == 255  # mask is where band == 255

            # Extract feature shapes and values from the array.
            rast_shapes = rasterio.features.shapes(band, transform=src.transform, mask=mask)
            geoms_geojson = [{'properties': {'val': val}, 'geometry': geom} for i, (geom, val) in
                             enumerate(rast_shapes)]
            gdf_poly = GeoDataFrame.from_features(list(geoms_geojson), crs=from_epsg(image_epsg))

            if config.get_debug_mode():
                temp_file_list += [os.path.join(TEMPDIR, '{}_{}BlockBoundary_{}.shp'
                                                .format(filename, len(temp_file_list) + 1,
                                                        image_epsg))]

                save_geopandas_tofile(gdf_poly, temp_file_list[-1], overwrite=True)

        del band, rast_shapes, geoms_geojson, mask

    else:
        desc_poly = VectorDescribe(polygon_shapefile)

        if 'POLY' not in desc_poly.geometry_type.upper():
            raise ValueError('Invalid input data : a polygon shapefile is required')

        gdf_poly = desc_poly.open_geo_dataframe()

        if groupby is not None and groupby not in gdf_poly.columns:
            raise ValueError('Groupby column {} does not exist'.format(groupby))

        if groupby is not None:
            # change null/nones to blank string
            gdf_poly[groupby].fillna('', inplace=True)

        # reproject shapefile to raster
        poly_epsg = desc_poly.crs.epsg_number
        if poly_epsg != image_epsg:
            # we have to reproject the vector points to match the raster
            print('WARNING: Projecting points from {} to {}'.format(desc_poly.crs.epsg, image_epsg))
            gdf_poly.to_crs(epsg=image_epsg, inplace=True)
            poly_epsg = image_epsg

        del desc_poly

        # If a column name is specified, dissolve by this and create multi-polygons
        # otherwise dissolve without a column name at treat all features as one
        if len(gdf_poly) > 1:
            if groupby is not None:
                gdf_poly = gdf_poly.dissolve(groupby, as_index=False).copy()
            else:
                gdf_poly = GeoDataFrame(geometry=[gdf_poly.unary_union])

    # Loop through polygns features ----------------------------------------------------------------
    output_files = []
    for index, feat in gdf_poly.iterrows():
        loop_time = time.time()
        step_time = time.time()
        status = '{} of {}'.format(index + 1, len(gdf_poly))

        if groupby is not None:
            feat_name = re.sub('[^0-9a-zA-Z]+', '-', str(feat[groupby])).strip('-')
            if feat_name == '':
                feat_name = 'No-Name'
        else:
            feat_name = ''

        # From shapes create a blockgrid mask. -----------------------------------------------------
        blockgrid_memfile = MemoryFile()
        meta_block = rasterio.open(image_file).meta.copy()
        blockgrid, new_blockmeta = convert_polygon_feature_to_raster(feat, pixel_size)
        meta_block.update(count=1, **new_blockmeta)

        with blockgrid_memfile.open(**meta_block) as dest:
            dest.write(blockgrid, 1)

        del blockgrid

        if config.get_debug_mode():
            LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                'Feature to block_grid', status, feat_name,
                dur=str(timedelta(seconds=time.time() - step_time))))

            temp_file_list += [os.path.join(TEMPDIR, '{}_{}BlockGrid_{}_{}.tif'.format(
                filename, len(temp_file_list) + 1, pixel_size_str, feat_name))]

            save_in_memory_raster_to_file(blockgrid_memfile, temp_file_list[-1])

        # Clip Dataset. ---------------------------------------------------------------------
        geoms = [mapping(feat.geometry)]

        with rasterio.open(image_file) as src:
            crop_data, crop_transform = rio_mask(src, geoms, crop=True)

            for iband in band_nums:  # range(1, src.count + 1):
                if crop_data[iband - 1].min() == crop_data[iband - 1].max():
                    raise GeometryError("There is no overlap between your polygon features "
                                        "and selected band(s)")

            meta = src.meta.copy()
            meta.update({"driver"   : "GTiff",
                         "height"   : int(crop_data.shape[1]),
                         "width"    : int(crop_data.shape[2]),
                         "transform": crop_transform})

            clip_memfile = MemoryFile()
            with clip_memfile.open(**meta) as dest:
                dest.write(crop_data)
                # reapply band tags to store the index name.
                for iband in range(1, src.count + 1):
                    # image statistics have changed so don't copy the tags
                    cleaned_tags = dict(
                        [(key, val) for key, val in six.iteritems(src.tags(iband)) if
                         not key.upper().startswith('STATISTIC')])
                    if len(cleaned_tags) > 0:
                        dest.update_tags(iband, **cleaned_tags)

        if config.get_debug_mode():
            LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                'Clipped Image to Feature', status, feat_name,
                dur=str(timedelta(seconds=time.time() - step_time))))

            temp_file_list += [os.path.join(TEMPDIR, '{}_{}crop_{}.tif'.format(
                filename, len(temp_file_list) + 1, feat_name))]

            save_in_memory_raster_to_file(clip_memfile, temp_file_list[-1])

        step_time = time.time()
        del crop_data, crop_transform, geoms

        # Resample to new pixel size using averaging -----------------------------------------------
        with clip_memfile.open() as src:
            meta = src.meta.copy()
            meta.update({'height'   : meta_block['height'],
                         'width'    : meta_block['width'],
                         'transform': meta_block['transform']})

            # Resample all bands and write to file
            resamp_memfile = MemoryFile()
            with resamp_memfile.open(**meta) as dest:
                for iband in range(1, src.count + 1):
                    reproject(source=rasterio.band(src, iband),
                              src_transform=src.transform,
                              destination=rasterio.band(dest, iband),
                              dst_transform=dest.transform,
                              resampling=Resampling.average)

                    # image statistics have changed so don't copy the tags
                    cleaned_tags = dict(
                        [(key, val) for key, val in
                         six.iteritems(src.tags(iband))
                         if not key.upper().startswith('STATISTIC')]
                    )
                    if len(cleaned_tags) > 0:
                        dest.update_tags(iband, **cleaned_tags)

        if config.get_debug_mode():
            LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                'Resampled to {}m'.format(pixel_size), status, feat_name,
                dur=str(timedelta(seconds=time.time() - step_time))))

            temp_file_list += [os.path.join(TEMPDIR, '{}_{}resampleAV_{}_{}.tif'.format(
                filename, len(temp_file_list) + 1, pixel_size_str, feat_name))]

            save_in_memory_raster_to_file(resamp_memfile, temp_file_list[-1])

        step_time = time.time()

        del clip_memfile

        # get a count of holes to see if fill is required.
        fillholes_memfile = MemoryFile()
        with resamp_memfile.open() as src:

            apply_nodata = int(src.nodata) if src.nodata.is_integer() else src.nodata

            with blockgrid_memfile.open() as src_bg:
                # find holes where blockgrid has data and resample is nodata - this can be used as
                # the mask in the fillnodata at a later stage if required.
                hole_count = np.logical_and((src_bg.read(1, masked=True).mask == 0),
                                            (src.read(1, masked=True).mask == 1)).sum()

            if hole_count == 0:
                with fillholes_memfile.open(**src.meta) as dest:
                    for iband in range(1, src.count + 1):
                        band = src.read(iband)

                        # reapply block grid mask
                        with blockgrid_memfile.open() as src_bg:
                            new_band = np.where((src_bg.read(1, masked=True).mask == 0),
                                                band, apply_nodata)

                        dest.write(new_band, iband)

                        # image statistics have changed so don't copy the tags
                        cleaned_tags = dict([
                            (key, val) for key, val in
                            six.iteritems(src.tags(iband))
                            if not key.upper().startswith('STATISTIC')
                        ])

                        if len(cleaned_tags) > 0:
                            dest.update_tags(iband, **cleaned_tags)

            # Fill holes ---------------------------------------------------------------------------
            if hole_count > 0:
                with fillholes_memfile.open(**src.meta) as dest:
                    for iband in range(1, src.count + 1):
                        # for rasterio 1.0.3+ the mask= options will interpolate values for all
                        # designated nodata pixels (marked by zeros in `mask`)
                        fill_nd = fillnodata(src.read(iband, masked=True),
                                             mask=None,  # Dont use mask option
                                             max_search_distance=100,  # default is 100
                                             smoothing_iterations=0)

                        # reapply block grid mask
                        with blockgrid_memfile.open() as src_bg:
                            filled = np.where((src_bg.read(1, masked=True).mask == 0),
                                              fill_nd, apply_nodata)

                        dest.write(filled, iband)
                        # image statistics have changed so don't copy the tags
                        cleaned_tags = dict([(key, val) for key, val in six.iteritems(src.tags(iband))
                                             if not key.upper().startswith('STATISTIC')])
                        if len(cleaned_tags) > 0:
                            dest.update_tags(iband, **cleaned_tags)

                        del fill_nd, filled

            if config.get_debug_mode():
                LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                    'Filled {} holes'.format(hole_count), status, feat_name,
                    dur=str(timedelta(seconds=time.time() - step_time))))

                step_time = time.time()
                temp_file_list += [os.path.join(TEMPDIR, '{}_{}filled_{}_{}.tif'.format(
                    filename, len(temp_file_list) + 1, pixel_size_str, feat_name))]

                save_in_memory_raster_to_file(fillholes_memfile, temp_file_list[-1])

        del resamp_memfile

        # Apply Focal Statistics to Smooth ---------------------------------------------------------
        with fillholes_memfile.open() as src:
            meta = src.meta.copy()

            if config.get_debug_mode():
                temp_file_list += [os.path.join(TEMPDIR, '{}_{}smoothed_{}_{}.tif'.format(
                    filename, len(temp_file_list) + 1, pixel_size_str, feat_name))]

            for i, iband in enumerate(band_nums, start=1):  # range(1, src.count + 1):
                smooth, _ = focal_statistics(src, iband, size=5, clip_to_mask=True)

                # change np.nan to a number
                smooth[np.isnan(smooth)] = apply_nodata
                meta.update(dtype=rasterio.dtypes.get_minimum_dtype(smooth))

                if config.get_debug_mode():
                    meta.update(count=len(band_nums))
                    with rasterio.open(temp_file_list[-1], 'w+', **meta) as dest:
                        dest.write(smooth, i)
                        # image statistics have changed so don't copy the tags
                        cleaned_tags = dict([(key, val) for key, val in src.tags(i).items()
                                             if not key.upper().startswith('STATISTIC')])
                        if len(cleaned_tags) > 0:
                            dest.update_tags(i, **cleaned_tags)

                if 'name' in src.tags(iband):
                    index_str = src.tags(iband)['name']
                else:
                    index_str = 'Band{}'.format(iband)

                if feat_name == '':
                    out_imagefile = '{}_{}'.format(index_str, pixel_size_str)
                else:
                    out_imagefile = '{}_{}_{}'.format(feat_name, index_str, pixel_size_str)

                out_imagefile = re.sub('-+', '-', out_imagefile).strip('_')
                out_imagefile = os.path.join(out_folder, out_imagefile + '.tif')

                meta.update(count=1)
                with rasterio.open(out_imagefile, 'w', tfw='YES', **meta) as dest, \
                        rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
                    dest.write(smooth, 1)

                    # image statistics have changed so don't copy the tags
                    cleaned_tags = dict([
                        (key, val) for key, val in
                        six.iteritems(src.tags(iband))
                        if not key.upper().startswith('STATISTIC')
                    ])
                    if len(cleaned_tags) > 0:
                        dest.update_tags(1, **cleaned_tags)

                output_files += [out_imagefile]

                del smooth

            if config.get_debug_mode():
                LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                    'Smoothed and saved to {} file(s)'.format(len(band_nums)),
                    status, feat_name, dur=str(timedelta(seconds=time.time() - step_time))))

                step_time = time.time()

            del fillholes_memfile, blockgrid_memfile

        LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
            'Created {} files for feature'.format(len(band_nums)),
            status, feat_name, dur=str(timedelta(seconds=time.time() - loop_time))))

        # del dest, src, src_bg

    return output_files


def calc_indices_for_block(image_file, pixel_size, band_map, out_folder, indices=[], image_epsg=0,
                           image_nodata=None, polygon_shapefile=None, groupby=None, out_epsg=0):
    """ Calculate indices for a multi band image then resample to a specified pixel size
      and block grid extent for each shapefile polygon.

      Use this tool to create single band images, one for each index and shapefile polygon
      combination. A group-by column may be used to dissolve multiple polygons belonging to an
      individual block.

      If a polygon shapefile is not specified, polygons will be created from the images' mask

      The polygon shapefile will be re-projected to match the image file

      A block grid will be created for each feature for the nominated pixel size and used as the
      base for analysis

      "image_epsg" and "image_nodata" can be used to set the coordinate system and image nodata
       values when they are not present within the image file.

      The output filename will consist of the selected feature and calculated index.

      The processing steps to achieve this are:
        - Reproject image to new coordinate system if required.
        - Calculate indices
        - Dissolve polygons optionally using the group-by column

        Loop through each polygon feature and......
          - Create Block Grid
          - Clip image to polygon
          - Resample and fit to block grid for a pixel size and using Average resampling technique
          - Identify holes and fill if necessary
          - smooth image using a 5x5 pixel moving average (focal_statistics)

    Args:
        image_file (str): the input image file
        pixel_size (int): the pixel size used for resampling
        band_map (pyprecag.bandops.BandMapping): A dictionary matching band numbers to
                                                band type (ie Red, Green, Blue etc.)
        out_folder (str): The output folder for the created images.
        indices (List[str]): The list of indices to calculate.
        image_epsg (int):  epsg number of the image to be used when missing in image
        image_nodata (int): nodata value of the image to be used when missing in image
        polygon_shapefile (str): a polygon shapefile used to cut up an image.
        groupby (str):  the column/field to use to group multiple features.
        out_epsg (int): The epsg number representing the coordinate system of the output images.

    Returns:
        List[str]: the list of created images.
    """

    if isinstance(polygon_shapefile, str) and polygon_shapefile.strip() == '':
        polygon_shapefile = None

    for ea_arg in [('image_file', image_file), ('out_folder', out_folder),
                   ('polygon_shapefile', polygon_shapefile)]:
        if ea_arg[0] != 'polygon_shapefile':
            if ea_arg[1] is None or ea_arg[1].strip() == '':
                raise ValueError('{} is required '.format(ea_arg[0]))

        if ea_arg[1] is not None and not os.path.exists(ea_arg[1]):
            raise IOError('{} does not exist-Got {}'.format(ea_arg[0], os.path.dirname(ea_arg[1])))

    for ea_arg in [('image_epsg', image_epsg), ('out_epsg', out_epsg)]:
        if not isinstance(ea_arg[1], six.integer_types):
            raise TypeError('{} must be a integer - Got {}'.format(*ea_arg))

    filename, ext = os.path.splitext(os.path.basename(image_file))

    if not os.path.basename(out_folder) == os.path.basename(image_file).replace('.', '_'):
        out_folder = os.path.join(out_folder, os.path.basename(image_file).replace('.', '_'))

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    with NamedTemporaryFile(prefix='{}_reproj_'.format(filename),
                            suffix='.tif', dir=TEMPDIR) as new_file:
        reproj_image = os.path.normpath(new_file.name)

    reproject_image(image_file, reproj_image, out_epsg, image_nodata=image_nodata,
                    image_epsg=image_epsg)

    with NamedTemporaryFile(prefix='{}_indices_'.format(filename),
                            suffix='.tif', dir=TEMPDIR) as new_file:
        indices_image = os.path.normpath(new_file.name)

    calculate_image_indices(reproj_image, band_map, indices_image, indices)

    if polygon_shapefile is None or polygon_shapefile.strip() == '':
        with NamedTemporaryFile(prefix='{}_BlockPolygons_'.format(filename),
                                suffix='.shp', dir=TEMPDIR) as new_file:
            polygon_shapefile = os.path.normpath(new_file.name)

        # create a polygon from the image. hopefully by now the nodata val is correct.
        with rasterio.open(reproj_image) as src:
            band = src.dataset_mask()
            mask = band == 255  # mask is where band == 255

            # Extract feature shapes and values from the array.
            rast_shapes = rasterio.features.shapes(band, transform=src.transform, mask=mask)
            geoms_geojson = [{'properties': {'val': val}, 'geometry': geom} for i, (geom, val) in
                             enumerate(rast_shapes)]
            gdf_poly = GeoDataFrame.from_features(list(geoms_geojson), crs=from_epsg(image_epsg))
            save_geopandas_tofile(gdf_poly, polygon_shapefile, overwrite=True)

        del band, rast_shapes, geoms_geojson, mask

    out_files = multi_block_bands_processing(indices_image, pixel_size, out_folder,
                                             polygon_shapefile=polygon_shapefile, groupby=groupby)

    if not config.get_debug_mode():
        if TEMPDIR in indices_image:
            os.remove(indices_image)

        if TEMPDIR in reproj_image:
            os.remove(reproj_image)

        if TEMPDIR in polygon_shapefile:
            os.remove(polygon_shapefile)

    return out_files


def resample_bands_to_block(image_file, pixel_size, out_folder, band_nums=[], image_epsg=0,
                            image_nodata=None, polygon_shapefile=None, groupby=None, out_epsg=0):
    """Derive multiple resampled image bands matching the specified pixel size and block grid extent
     for each shapefile polygon.

    Use this tool create individual band images for each polygon within a shapefile. A group-by
    column may be used to dissolve multiple polygons belonging to an individual block. The fitting
    of rasters to a base Block (grid) ensures for easier, more accurate multi-layered analysis
    required by in Precision Agriculture.

    The processing steps to achieve this are:
        - Reproject image to nominated coordinate system
        - Dissolve polygons optionally using the groupby column

        Loop through each polygon feature and......
          - Create Block Grid
          - Clip image to polygon
          - Resample and fit to block grid for a pixel size and using Average resampling technique
          - Identify holes and fill if necessary
          - smooth image using a 5x5 pixel moving average (focal_statistics)

    If a polygon shapefile is not specified, polygons will be created from the images' mask

    The polygon shapefile will be re-projected to match the image file

    A block grid will be created for each feature for the nominated pixel size and used as the
    base for analysis

    "image_epsg" and "image_nodata" can be used to set the coordinate system and image nodata values
     when they are not present within the image file.

    The output filename will consist of the selected band number or the value of a band's rasterio
    custom name tag

    Args:
        image_file (str): An input image
        pixel_size (int, float): The desired output pixel size in metres.
        out_folder (str): The output folder for the created images.
        band_nums  (List[int]): a list of band numbers to process. If empty all bands will be used
        image_epsg (int):  epsg number of the image to be used when missing in image
        image_nodata (int): nodata value of the image to be used when missing in image
        polygon_shapefile (str): a polygon shapefile used to cut up an image.
        groupby (str):  the column/field to use to group multiple features.
        out_epsg (int): The epsg number representing the output coordinate system.

    Returns:
        List[str]: a list of created files.
    """

    filename, ext = os.path.splitext(os.path.basename(image_file))

    for ea_arg in [('image_file', image_file), ('out_folder', out_folder)]:
        if ea_arg[1] is None or ea_arg[1].strip() == '':
            raise ValueError('{} is required '.format(ea_arg[0]))

        if ea_arg[1] is not None and not os.path.exists(ea_arg[1]):
            raise IOError('{} does not exist-Got {}'.format(ea_arg[0], os.path.dirname(ea_arg[1])))

    if not os.path.basename(out_folder) == os.path.basename(image_file).replace('.', '_'):
        out_folder = os.path.join(out_folder, os.path.basename(image_file).replace('.', '_'))

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    with NamedTemporaryFile(prefix='{}_reproj_'.format(filename),
                            suffix='.tif', dir=TEMPDIR) as new_file:
        reproj_image = os.path.normpath(new_file.name)

    reproject_image(image_file, reproj_image, out_epsg, image_nodata=image_nodata,
                    image_epsg=image_epsg)

    out_files = multi_block_bands_processing(reproj_image, pixel_size, out_folder, groupby=groupby,
                                             band_nums=band_nums,
                                             polygon_shapefile=polygon_shapefile)

    if TEMPDIR in reproj_image and not config.get_debug_mode():
        os.remove(reproj_image)

    return out_files


def kmeans_clustering(raster_files, output_tif, n_clusters=3, max_iterations=500):
    """Create zones with k-means clustering from multiple raster files as described in


    The input raster files should all:
        - have the same pixel size
        - be in the same coordinate system
        - should overlap
    Only the first band of each raster will be used.

    The output TIFF image extent will be the minimum overlapping extent of the input images.
    Each image will be resampled to a fixed coordinate to ensure pixels between images align.

    Args:
        raster_files (dict[str]): The dictionary ({filename:alias}) of input raster files
                                  and alias names.
        output_tif (str):   The output TIFF file
        n_clusters (int):  The number of clusters/zones to create.
        max_iterations (int): Maximum number of iterations for k-means algorithm in a
                              single run.

    Returns:
        pandas.core.frame.DataFrame: A dataframe containing cluster statistics for each image.
    """

    if not isinstance(n_clusters, six.integer_types):
        raise TypeError('Size must be an Integer.')

    if not isinstance(max_iterations, six.integer_types):
        raise TypeError('Size must be an Integer.')

    if not isinstance(raster_files, (list, dict)):
        raise TypeError('Invalid Type: raster_files should be a dictionary of filename and alias')

    if isinstance(raster_files, list):
        warnings.warn('raster_files has now been converted to dictionary pairs of filenames and unique alias. \n '
                      'Converting now', DeprecationWarning)
        raster_files = {ea: os.path.basename(os.path.splitext(ea)[0]) for ea in raster_files}
    else:
        if len(set(raster_files.values())) == 1:  # all aliases are the same.
            warnings.warn('raster_files have no aliases, using filename')
            raster_files = {ea: os.path.basename(os.path.splitext(ea)[0]) for ea in raster_files}
        elif len(set(raster_files.values())) != len(raster_files.values()):
            raise ValueError('Raster Aliases are not unique')

    not_exists = [my_file for my_file in raster_files.keys() if not os.path.exists(my_file)]
    if len(not_exists) > 0:
        raise IOError('raster_files: {} raster file(s) do '
                      'not exist\n\t({})'.format(len(not_exists), '\n\t'.join(not_exists)))

    temp_file_list = []
    start_time = time.time()

    images_combined, band_list = stack_and_clip_rasters(list(raster_files.keys()),
                                                        output_tif='1_kmeans_combined.tif')
    temp_file_list.append(images_combined)

    # Make sure ALL masks are saved inside the TIFF file not as a sidecar.
    gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', 'YES')

    step_time = time.time()

    # Extract data ready for k-means ------------------------------------------------
    with rasterio.open(images_combined) as src:
        template_band = src.read(1, masked=True)

        stack_meta = src.meta.copy()

        bands = src.read(masked=True)
        data_stack = []
        for ea_band in bands:
            new_band = np.ma.masked_values(ea_band, src.nodata)
            data_stack.append(np.ma.compressed(ea_band))

        data_stack = np.array(data_stack)

    # Now we need to whiten and/or normalise the data prior to clustering
    # This ensures each variable carries the same weight despite different input data ranges
    norm_stack = np.array([stats.zscore(ea) for ea in data_stack])

    if config.get_debug_mode():
        # write mask to file
        temp_file_list += [os.path.join(
            TEMPDIR, '{}_kmeans_normalised_image.tif'.format(len(temp_file_list) + 1))]

        with rasterio.open(temp_file_list[-1], 'w', **stack_meta) as tmp_dst:

            for i, ea in enumerate(norm_stack, start=1):
                # recreate the full array dimensions
                new_band = np.empty_like(template_band)
                np.place(new_band, ~new_band.mask, ea)

                # reapply the mask to remove values assigned to nodata pixels
                new_band = np.ma.masked_values(new_band, tmp_dst.nodata)
                tmp_dst.update_tags(i, name=band_list[i - 1])

                tmp_dst.write(new_band, i)
            tmp_dst.descriptions = band_list

        LOGGER.info('{:<30} {:<15} {dur}'.format('Images Normalised', '',
                                                 dur=str(timedelta(seconds=time.time() - step_time))))
        step_time = time.time()

    # Run the k-means clustering -------------------------------------------------------------------

    _, code = kmeans2(data=norm_stack.T, k=n_clusters, iter=max_iterations, minit='points')

    # adjust the class values to start from 1 so 0 can be used as nodata

    # Write to file
    cluster_data = np.empty_like(template_band, dtype='int')
    np.place(cluster_data, ~cluster_data.mask, code.T)

    # set nodata to 0
    cluster_data = np.ma.masked_values(cluster_data, -1)
    cluster_dtype = rasterio.dtypes.get_minimum_dtype([-1] + cluster_data.data)

    cluster_meta = stack_meta.copy()
    cluster_meta.update({'count': 1, 'nodata': -1, 'dtype': cluster_dtype})

    # create a raster in memory
    clust_memfile = MemoryFile()
    # open the in-memory raster for writing
    with clust_memfile.open(**cluster_meta) as dest:
        dest.write(cluster_data.astype(cluster_dtype), 1)

    LOGGER.info('{:<30} {:<15} {dur}'.format('Clustering complete', '',
                                             dur=str(timedelta(seconds=time.time() - step_time))))
    step_time = time.time()

    # create summary statistics --------------------------------------------------------------------

    # create a table to store statistics
    results_df = pd.DataFrame(columns=['zone'])

    with clust_memfile.open() as src_clust, \
            rasterio.open(images_combined) as src_img:

        # get the list of clusters to process
        clust_list = list(np.unique(src_clust.read()))
        clust_list.remove(src_clust.nodata)

        bands = src_img.read(masked=True)

        for ea_clust in clust_list:
            # get pixels for each cluster
            clust_mask = np.where((src_clust.read(1, masked=True) == ea_clust), 1, 0)

            new_row = pd.DataFrame([ea_clust], columns=['zone'])
            img_meta = src_img.meta.copy()

            with MemoryFile() as memfile:
                # apply cluster mask to all bands
                with memfile.open(**img_meta) as tmp_dst:
                    tmp_dst.write_mask(clust_mask.astype(img_meta['dtype']))
                    tmp_dst.write(bands)

                with memfile.open() as tmp_src:
                    for i, ea_band in enumerate(tmp_src.read(masked=True)):
                        alias = next(val for key, val in raster_files.items() if src_img.descriptions[i] == Path(key).stem)

                        if not alias:
                            alias = src_img.descriptions[i]

                        new_row[alias + ', mean'] = np.nanmean(ea_band)
                        new_row[alias + ', std'] = np.nanstd(ea_band)

                        # add a blank column inorder for populating later
                        new_row[alias + ', vesper'] = np.nan

            results_df = pd.concat([results_df, new_row])

        # reorder zone numbering based on mean value of all inputs
        mean_cols = [col for col in results_df.columns if ', mean' in col]

        results_df['zone_mean'] = results_df[mean_cols].mean(axis=1)
        results_df['new_zone'] = results_df['zone_mean'].rank(ascending=False).astype(int)
        n_dict = results_df.set_index('zone').to_dict()['new_zone']

        # remap it in the statistics
        results_df['zone'] = results_df['new_zone']
        results_df.drop(['zone_mean', 'new_zone'], axis=1, inplace=True)
        results_df.sort_values(by='zone', inplace=True)

        # add nodata to 0 to remap
        n_dict[cluster_data.fill_value] = 0

        # remap to new values
        reclass = cluster_data.copy()
        # set the nodata value
        reclass.fill_value = 0
        reclass.filled(0)

        # remap the array
        for k, v in n_dict.items():
            reclass.data[cluster_data.data == k] = v

        # update the metadata
        cluster_dtype = rasterio.dtypes.get_minimum_dtype([reclass.fill_value] + reclass.data)
        cluster_meta.update({'nodata': reclass.fill_value,
                             'dtype' : cluster_dtype})

        # and write to file
        with rasterio.open(output_tif, 'w', **cluster_meta) as dst:
            dst.write(reclass.astype(cluster_dtype), 1)

        # loop through and add 95Conf level & MeanPredSE
        new_row = pd.DataFrame(['95ConfLevel', 'MedianPredSE'], columns=['zone'])

        col_order = list(results_df.columns)
        col_order.insert(0, col_order.pop(col_order.index('zone')))

        with rasterio.open(images_combined) as src_img:
            for iband in range(1, src_img.count + 1):
                for tag in ['PAT_MedianPredSE', 'PAT_95ConfLevel']:
                    if tag in src_img.tags(iband):
                        tag_col = tag.replace('PAT_', '')
                        # get the alias
                        alias = next(val for key, val in raster_files.items() if src_img.tags(iband)['name'] in key)
                        if not alias:
                            alias = src_img.tags(iband)['name']

                        new_row.loc[new_row['zone'] == tag_col, '{}, vesper'.format(alias)] = src_img.tags(iband)[tag]

        results_df = pd.concat([results_df, new_row])

        results_df.dropna(axis='columns', how='all', inplace=True)

        # write to csv without ', ' to assist when loading into ESRI
        col_names = results_df.columns.str.replace(', ', '_').values
        results_df.to_csv(output_tif.replace('.tif', '_statistics.csv'),
                          header=col_names, index=False)

        # format the table and print to log.
        results_df_copy = results_df.copy()
        col_names = results_df_copy.columns.str.split(', ', expand=True).values

        # replace column name NaNs to '....'. '....' as multiple spaces aren't allowed
        results_df_copy.columns = pd.MultiIndex.from_tuples(
            [('.......', x[0]) if pd.isnull(x[1]) else x for x in col_names])

        LOGGER.info('\nCluster Statistics:\n{}\n'.format(
            results_df_copy.to_string(justify='center', index=False)))

        del results_df_copy

    LOGGER.info('{:<30} {:<15} {dur:<15} {}'.format(
        'Saved Stats CSV', '', output_tif.replace('.tif', '_statistics.csv'),
        dur=str(timedelta(seconds=time.time() - step_time))
    ))

    # clean up of intermediate files
    if len(temp_file_list) > 0 and not config.get_debug_mode():
        for ea in temp_file_list:
            for f in glob.glob(os.path.splitext(ea)[0] + '*'):
                os.remove(f)

    gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', None)

    LOGGER.info('{:<30}  {:<15} {dur:<15} {}'.format(
        'K-Means Clustering Completed', '', '{} zones for {} rasters'.format(
            n_clusters, len(raster_files)),
        dur=str(timedelta(seconds=time.time() - start_time))
    ))

    return results_df


def create_points_along_line(lines_geodataframe, lines_crs, distance_between_points,
                             offset_distance, out_epsg=0, out_points_shapefile=None,
                             out_lines_shapefile=None):
    """Add points along a line using a specified distance and create left/right parallel points
        offset by a distance.

        If the lines are in a geographic coordinate system they will be re-projected to a projected
        coordinate system.

        All touching lines will be treated as one.
        MultiPart geometry will be converted to single part geometry.
        The first and last points will be offset from start/end of the line evenly.
        Attributes from the input lines will be lost.

        line_crs is used to ensure that the correct wkt definition is maintained when using
        geopandas.

        Args:
            lines_geodataframe (geopandas.geodataframe.GeoDataFrame): A Geopandas dataframe
                             containing Lines
            lines_crs (pyprecag.crs.crs): The detailed coordinate system

            distance_between_points (int): The separation distance between points.

            offset_distance (int): The distance between the Strip point and parallel point.

            out_epsg (int): Optionally specify the epsg number for the output coordinate system.
                            This should be a project coordinate system

            out_points_shapefile (str): Optionally specify shapefile path and filename used to
                            save the points to. If a path is not supplied, it will save the file
                            to TEMPDIR  by default.

            out_lines_shapefile (str): Optionally specify shapefile path and filename used to save
                            the lines to. If a path is not supplied, it will save the file to
                            TEMPDIR  by default
        Returns:
             geopandas.geodataframe.GeoDataFrame: The geodataframe containing the created points.
             pyprecag.crs.crs: The coordinate system of both the points and lines geodataframe.
             geopandas.geodataframe.GeoDataFrame: The geodataframe containing the created lines.
        """

    if not isinstance(lines_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : inputGeodataFrame')

    if not any("LINE" in g.upper() for g in lines_geodataframe.geom_type.unique()):
        raise GeometryError('Invalid input data : A lines geopandas dataframe is required')

    for argCheck in [('offset_distance', offset_distance),
                     ('distance_between_points', distance_between_points)]:
        if not isinstance(argCheck[1], number_types):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not isinstance(lines_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if not isinstance(out_epsg, six.integer_types):
        raise TypeError('out_epsg must be a integer - Got {}'.format(*out_epsg))

    if out_points_shapefile is not None and os.path.isabs(out_points_shapefile):
        if not os.path.exists(os.path.dirname(out_points_shapefile)):
            raise IOError(
                'Output directory {} does not exist'.format(os.path.dirname(out_points_shapefile)))

    if out_lines_shapefile is not None and os.path.isabs(out_lines_shapefile):
        if not os.path.exists(os.path.dirname(out_lines_shapefile)):
            raise IOError(
                'Output directory {} does not exist'.format(os.path.dirname(out_lines_shapefile)))

    if out_points_shapefile is not None and out_lines_shapefile is not None:
        if out_points_shapefile == out_lines_shapefile:
            raise IOError('Output points and lines shapefile names are identical')

    # create temp files template
    with NamedTemporaryFile(prefix="strip_treatment_", suffix='.shp', dir=TEMPDIR) as new_file:
        temp_filename = new_file.name

    start_time = time.time()
    step_time = time.time()

    # input needs to be a projected coordinate system to work with metric distances
    if lines_geodataframe.crs.is_geographic:
        if out_epsg <= 0:
            xmin, ymin, _, _ = lines_geodataframe.total_bounds
            points_crs = pyprecag_crs.getProjectedCRSForXY(xmin, ymin, lines_crs.epsg_number)
    else:
        points_crs = lines_crs

    # project if required.
    if lines_crs.epsg_number != points_crs.epsg_number:
        lines_geodataframe.to_crs(epsg=points_crs.epsg_number, inplace=True)

        if config.get_debug_mode():
            LOGGER.info('{:<30}   {:<15} {dur}'.format('Reproject lines To epsg {}'.format(points_crs.epsg_number), '',
                                                       dur=str(timedelta(seconds=time.time() - step_time))))

    step_time = time.time()

    # merge touching lines into single parts
    gdf_lines = GeoDataFrame(geometry=[lines_geodataframe.geometry.agg(lambda x: linemerge(x.tolist()))],
                             crs=lines_geodataframe.crs)
    gdf_lines = gdf_lines.explode(index_parts=False).reset_index(drop=False, names='FID')

    # Add TrialID  Side, Length and startoffset
    if 'TrialID' not in gdf_lines.columns:
        gdf_lines.insert(0, 'TrialID', gdf_lines.index)
    else:
        gdf_lines['TrialID'] = gdf_lines.index

    # Add Side, Length and startoffset attributes
    gdf_lines['name_s'] = 'Strip'
    gdf_lines['bearing'] = gdf_lines.geometry.apply(line_bearing)

    # Calculate the length of each line
    gdf_lines['length'] = gdf_lines['geometry'].length
    
    ## Add Offsets
    gdf_lines['line_r'] = gdf_lines.apply(lambda l: l.geometry.offset_curve(-offset_distance, quad_segs=16, join_style='mitre'), axis=1)
    gdf_lines['line_l'] = gdf_lines.apply(lambda l: l.geometry.offset_curve(offset_distance, quad_segs=16, join_style='mitre'), axis=1)
    gdf_lines['name_r'] = gdf_lines['line_r'].apply(lambda x: deg_to_8_compass_pts(line_bearing(x) + 90) + ' Strip')
    gdf_lines['name_l'] = gdf_lines['line_l'].apply(lambda x: deg_to_8_compass_pts(line_bearing(x) - 90) + ' Strip')

    striplines_df = pd.wide_to_long(gdf_lines.drop(columns=['length','bearing','FID']).rename(columns={'geometry': 'line_s'}),
                                    stubnames=['line', 'name'], sep='_',
                                    i=['TrialID'], j='Strip', suffix=r'\w+')
    striplines_df = striplines_df.reset_index().rename(columns={'name': 'Strip_Name', 'line': 'geometry'})

    striplines_gdf = GeoDataFrame(striplines_df.drop(columns='Strip'), crs=gdf_lines.crs)
    del striplines_df
    def split_line_by_distance(line, split_distance, center_on_line=True):
        start_offset = 0
        if center_on_line:
            start_offset = (line.length % split_distance) / 2

        start_pts = np.arange(0+start_offset, line.length, split_distance)
        end_pts = np.roll(start_pts, -1)

        seg_lines = [substring(force_2d(line), start_dist=s, end_dist=e, normalized=False)
                     for s, e in zip(start_pts[:-1], end_pts[:-1])]

        return MultiLineString(seg_lines)

    def make_paired_segments(line, paired_line):
        start = paired_line.project(Point(line.coords[0]))
        end = paired_line.project(Point(line.coords[-1]))
        return substring(paired_line, start_dist=start, end_dist=end, normalized=False)

    gdf_lines['seg_s'] = gdf_lines.geometry.apply(split_line_by_distance,
                                                  split_distance=distance_between_points,
                                                  center_on_line=True)

    exp_df = gdf_lines.rename_geometry('line_s').set_geometry('seg_s').explode(index_parts=False)
    exp_df.reset_index(drop=True, inplace=True)

    exp_df.insert(1, 'PointID', exp_df.groupby('TrialID').cumcount())
    exp_df['length'] = exp_df.set_geometry('seg_s').length
    exp_df.insert(3, 'DistOnLine', exp_df.groupby(['TrialID'])['length'].cumsum())

    exp_df['seg_r'] = exp_df.apply(lambda l: make_paired_segments(l['seg_s'], l['line_r']), axis=1)
    exp_df['seg_l'] = exp_df.apply(lambda l: make_paired_segments(l['seg_s'], l['line_l']), axis=1)

    seg_cols = exp_df.filter(like='seg_').columns.to_list()
    name_cols = exp_df.filter(like='name_').columns.to_list()

    segs_gdf = pd.wide_to_long(exp_df[['TrialID', 'PointID', 'DistOnLine']+ seg_cols + name_cols],
                                   stubnames=['seg', 'name'], sep='_',
                                   i=['TrialID', 'PointID', 'DistOnLine'], j='Strip', suffix=r'\w+')
    segs_gdf = segs_gdf.reset_index().rename(columns={'name': 'Strip_Name', 'seg': 'geom_seglines'})

    segs_gdf['geometry'] = segs_gdf.set_geometry('geom_seglines').centroid
    segs_gdf = GeoDataFrame(segs_gdf,crs=gdf_lines.crs)

    if config.get_debug_mode():
        segs_gdf['length'] = segs_gdf['geom_seglines'].length
        save_geopandas_tofile(segs_gdf.set_geometry('geom_seglines', drop=True), temp_filename.replace('.shp', '_seglines.shp'), overwrite=True)

    segs_gdf.drop(columns=['Strip','geom_seglines'],inplace=True)

    if out_lines_shapefile is not None or config.get_debug_mode():
        if out_lines_shapefile is None:
            save_geopandas_tofile(striplines_gdf, temp_filename.replace('.shp', '_lines.shp'), overwrite=True)
        else:
            save_geopandas_tofile(striplines_gdf, out_lines_shapefile, overwrite=True)

    if out_points_shapefile is not None or config.get_debug_mode():
        if out_points_shapefile is None:
            save_geopandas_tofile(segs_gdf, temp_filename.replace('.shp', '_points.shp'), overwrite=True)
        else:
            save_geopandas_tofile(segs_gdf, out_points_shapefile, overwrite=True)

    if config.get_debug_mode():
        LOGGER.info('{:<30} {:>15} {dur}'.format('Create Points Along Line Completed', '',
                                                 dur=str(timedelta(seconds=time.time() - start_time))))

    return segs_gdf, points_crs, striplines_gdf


def ttest_analysis(points_geodataframe, points_crs, values_raster, out_folder,
                   zone_raster='', control_raster='', size=5, create_graph=False):
    """Run a moving window t-test analysis for a strip trial as described in Lawes
    and Bramley (2012).

    Format of the points must be from the create_points_along_line tools.
    All input rasters must be of the same coordinate system and pixel size and overlap with the
    points.

    Output statistics include:
        controls_mean  - row by row mean of the control columns
        treat_diff -   row by row difference between the treatment and controls_mean columns
        av_treat_diff - calculate mean of values using a moving window using the treat_diff column
        p_value - calculate  p_value using a moving window using treatment and controls_mean columns
        RI  - Response Index using the treatment and controls_mean columns

    Output Files include:
        For each line and strip combination :
            - png Map showing orientation of the line (start and ends)
            - png set of graphs
            - CSV file of derived statistics.
    Reference:
        Lawes RA, Bramley RGV. 2012. A Simple Method for the Analysis of On-Farm Strip Trials.
         Agronomy Journal 104, 371-377.

    Args:
        points_geodataframe (geopandas.geodataframe.GeoDataFrame): points derived using
                create_points_along_line
        points_crs (pyprecag.crs.crs):the coordinate system for the points.
        values_raster (str): a kriged raster containing the treatment  values
        out_folder (str):  folder for output files.
        zone_raster (str): a raster containing zones.
        control_raster (str): a kriged raster of
        size (int): the size used to calculate the moving window statistics.
        create_graph (bool):

    Returns:
        pandas.core.frame.DataFrame: dataframe containing output statistics
    """
    warnings.warn("Moving window T-test will be modified to include Dutilleul's correction. For now the reported \
                   p value should not be regarded as a true indication of significance", DeprecationWarning)

    if not isinstance(points_geodataframe, GeoDataFrame):
        raise TypeError('Invalid input data : inputGeodataFrame')

    if not any("POINT" in g.upper() for g in points_geodataframe.geom_type.unique()):
        raise GeometryError('Invalid input data : a points geopandas dataframe is required')

    if not isinstance(points_crs, pyprecag_crs.crs):
        raise TypeError('Crs must be an instance of pyprecag.crs.crs')

    if out_folder is None or out_folder == '':
        raise ValueError('Please specify an output folder')

    if not os.path.exists(out_folder):
        raise IOError('Output directory {} does not exist'.format(out_folder))

    if not isinstance(size, int):  # or size % 2 == 0   - odd numbers only
        raise TypeError("Size {} should be an integer.".format(size))

    if zone_raster is None:
        zone_raster = ''

    if control_raster is None:
        control_raster = ''

    not_exists = [ea for ea in [values_raster, control_raster, zone_raster] if ea and
                  not os.path.exists(ea)]

    if len(not_exists) > 0:
        raise IOError('raster files: {} raster file(s) '
                      'do not exist\n\t({})'.format(len(not_exists), '\n\t'.join(not_exists)))
    del not_exists

    missing = [ea for ea in ['TrialID', 'PointID', 'Strip_Name', 'DistOnLine']
               if ea and ea not in points_geodataframe.columns]

    if len(missing) > 0:
        raise ValueError("columns not found - {}. Please use the output from the "
                         "create_points_along_line. (PAT's create strip trial points)"
                         " ".format(len(missing), ','.join(missing)))

    del missing

    # get list of valid rasters
    raster_files = [ea for ea in [values_raster, control_raster, zone_raster] if ea]

    # make sure rasters overlap  ------------------------------------
    try:
        for i, ea_raster in enumerate(raster_files, start=1):
            with rasterio.open(ea_raster) as src:

                band1 = src.read(1, masked=True)

                # get minimum data extent
                data_window = get_data_window(band1)

                # find the intersection of all the windows
                if i == 1:
                    min_window = data_window
                else:
                    # create a new window using the last coordinates based on this image in case
                    # extents/pixel origins are different
                    min_img_window = from_bounds(*min_bbox, transform=src.transform)

                    # find the intersection of the windows.
                    min_window = intersection(min_img_window, data_window).round_lengths()

                # convert the window co coordinates
                min_bbox = src.window_bounds(min_window)

                del min_window

    except rasterio.errors.WindowError as e:
        # reword 'windows do not intersect' error message
        if not e.args:
            e.args = ('',)
        e.args = ("Rasters do not overlap",)
        raise  # re-raise current exception

    transect_column = 'Strip_Name'
    x_axis_column = 'DistOnLine'

    line_count = len(points_geodataframe['TrialID'].unique())

    # Extract raster values for points ------------------------------------------------------------
    gdf_points, points_crs = extract_pixel_statistics_for_points(points_geodataframe, points_crs,
                                                                 raster_files, 'extract_pixels.csv',
                                                                 size_list=[1])

    gdf_points.index.name = 'rowID'
    column_names = {}

    # find the columns relating to the rasters
    if values_raster != '':
        rastname = os.path.splitext(os.path.basename(values_raster))[0]
        column_names['Value'] = next(eaString for eaString in gdf_points.columns
                                     if rastname in eaString)

    if control_raster != '':
        rastname = os.path.splitext(os.path.basename(control_raster))[0]
        column_names['Control'] = next(eaString for eaString in gdf_points.columns
                                       if rastname in eaString)

    if zone_raster != '':
        rastname = os.path.splitext(os.path.basename(zone_raster))[0]
        column_names['Zone'] = next(eaString for eaString in gdf_points.columns
                                    if rastname in eaString)

    else:
        # Add a zone col with a constant value
        gdf_points['Zone'] = 1
        column_names['Zone'] = 'Zone'

    # rename columns from raster name to the dictionary key to be more generic
    gdf_points.rename(columns=dict([[v, k] for k, v in six.iteritems(column_names)]), inplace=True)

    # group by trialid and find those where all are none ie no overlap with rasters
    dfnulls = gdf_points.set_index('TrialID')['Value'].isnull().groupby('TrialID').all()

    # get only the trial id's with data
    gdf_points = gdf_points[~gdf_points['TrialID'].isin(dfnulls[dfnulls].index.values.tolist())]
    del dfnulls

    # create a unique id for each line/point
    gdf_points.insert(1, 'TrialPtID', gdf_points.apply(lambda x: "'{}-{}'".format(x['TrialID'], x['PointID']), axis=1))

    # Apply bearing to all transects
    for group_name, gdf_group in gdf_points.groupby(['TrialID', transect_column]):
        gdf_group = gdf_group.copy()
        gdf_group['geom2'] = gdf_group['geometry'].shift()
        gdf_group['label_angle'] = gdf_group.dropna(subset=['geom2'], axis=0).apply(
            lambda p: text_rotation(p['geometry'], p['geom2']), axis=1)

        # update the main dataframe - needs matching indexes
        gdf_points.loc[gdf_points.index.isin(gdf_group.index), 'label_angle'] = \
            gdf_group[['label_angle']]

    del group_name, gdf_group

    offset_names = gdf_points['Strip_Name'].unique().tolist()
    strip_name = offset_names.pop(offset_names.index('Strip'))

    # pivot/transpose the table.
    df_pivot = pd.DataFrame(gdf_points).pivot(index='TrialPtID', columns=transect_column,
                                              values=column_names.keys())

    # clean up un required columns
    if control_raster != '':  # only keep offset strips
        dropcols = [('Control', strip_name)]
        dropcols += [('Value', ea) for ea in offset_names]
        df_pivot.drop(dropcols, axis=1, inplace=True)

    # change the order before joining levels together
    df_pivot.columns = df_pivot.columns.swaplevel(0, 1)

    # Fix column headers
    df_pivot.columns = df_pivot.columns.map(' '.join).str.strip()

    ''' Ensure zone columns are integer. cant be integer and contain NAN's so do it later.
     TODO: In v0.24, you can now do df = df.astype(pd.Int32Dtype())
     http://pandas.pydata.org/pandas-docs/stable/user_guide/integer_na.html
     df_pivot[column_names['Zone']] = df_pivot[column_names['Zone']].astype(pd.Int8Dtype())'''

    # Prepare to join tables --------------------------------------------------------------------
    # drop geometry etc. and create flat table.
    df_table = pd.DataFrame(gdf_points.drop(columns='geometry', inplace=False))

    # Keep only minimum of columns
    dropcols = [ea for ea in df_table.columns if
                ea not in ['FID', 'TrialPtID', 'TrialID', 'PointID', 'DistOnLine']]

    df_table.drop(dropcols, axis=1, inplace=True)

    # remove duplicate points only keeping first - ie the Strip
    df_table.drop_duplicates(subset=['TrialPtID'], keep='first', inplace=True)

    # add pivot columns back to original table --------------------------------------------------
    df_table = pd.merge(df_table, df_pivot, on='TrialPtID')

    pivot_columns = df_pivot.columns
    del df_pivot

    import matplotlib as mpl
    # suppress ImportError: No module named _tkinter, please install the python-tk package
    mpl.use('Agg')

    import matplotlib.pyplot as plt
    # Suppress interactive plotting
    plt.ioff()

    from matplotlib import gridspec
    from collections import OrderedDict
    import matplotlib.patheffects as pe

    for iline, (line_id, gdf_strip) in enumerate(gdf_points.groupby(by='TrialID'), start=1):
        status = '{} of {}'.format(iline, line_count)
        loop_time = time.time()
        df_subtable = df_table[df_table['TrialID'] == line_id].copy()

        offset_names = gdf_strip['Strip_Name'].unique().tolist()
        strip_name = offset_names.pop(offset_names.index('Strip'))

        # determine which columns to use
        keepcols = [col for col in pivot_columns if len(col.split(' ')) == 2]
        keepcols += [col for x in offset_names for col in pivot_columns if col.startswith(x)]
        dropcols = [ea for ea in pivot_columns if ea not in keepcols]

        column_names = defaultdict(list)
        for ea in keepcols:
            parts = ea.split(' ')
            if len(parts) == 2:
                column_names[parts[-1]].append(ea)
            elif parts[-1] != 'Zone':
                column_names['Control'].append(ea)
            else:
                column_names[parts[-1]].append(ea)

        df_subtable.drop(dropcols, axis=1, inplace=True)

        # Convert zones to integer - there should be no NAN's
        df_subtable[column_names['Zone']] = df_subtable[column_names['Zone']].astype(int)

        # run loop for the three scenarios. 1) both sides (N+S), 2) North, 3) South
        scenarios = [offset_names] + offset_names
        for iscenario, scenario in enumerate(scenarios, start=1):

            if not isinstance(scenario, list):
                scenario = [scenario]

            control_col = [col for x in scenario for col in column_names['Control']
                           if col.startswith(x)]
            valid_zone_col = [col for x in scenario + ['Strip'] for col in column_names['Zone']
                              if col.startswith(x)]

            dropcols = [col for col in column_names['Control']
                        if col not in control_col + valid_zone_col + column_names['Value']]
            dropcols += [col for col in column_names['Zone']
                         if col not in control_col + valid_zone_col + column_names['Value']]

            offset_subst = '-'.join([ea.split(' ')[0] for ea in control_col])
            file_path_noext = os.path.join(out_folder, "Trial-{}_{}-strip".format(line_id,
                                                                                  offset_subst))

            df_statstable, control_mean = calculate_strip_stats(df_subtable.drop(columns=dropcols),
                                                                column_names['Value'][0],
                                                                control_col, size=size)

            df_statstable.drop(columns=['FID', 'TrialPtID'], axis=1).to_csv(
                file_path_noext + '.csv', index=False)

            if config.get_debug_mode():
                LOGGER.info('{:<30}\t{:>10}   {dur:<15} {}'.format(
                    'Saved CSV', '', file_path_noext + '.csv',
                    dur=str(timedelta(seconds=time.time() - loop_time))))

            ''' Create the map & graphs ---------------------------------------------------------
            https://stackoverflow.com/a/39694347
            https://stackoverflow.com/a/14887119'''

            markers = ['.', '^', '+', 'x', '*', 'o', ',', 'v', '<', '>', 's', 'd']
            colour_list = ['r', 'g', 'b', 'y']
            zone_column = 'Strip Zone'

            # add plotting parameters to table -----------------------------------------------
            # assign an id to the column
            df_statstable['zone_UID'] = df_statstable.groupby(zone_column).grouper.group_info[0]

            # assign a marker based on the index from the list
            df_statstable["zone_marker"] = df_statstable['zone_UID'].apply(lambda x: markers[x])

            # drop nan's introduced by the rolling window
            df_statstable = df_statstable.dropna(subset=['av_treat_dif'], axis=0).copy()

            # only get those points geometry with rolling window results.
            gdf_map = pd.merge(gdf_strip, df_statstable[['RI']],
                               on='TrialPtID').set_index('FID', drop=False)  # ,'label_angle'

            # get y axis limits
            min_ylimit = df_statstable['av_treat_dif'].min(skipna=True)
            max_ylimit = df_statstable['av_treat_dif'].max(skipna=True)

            use_pvalue = 0.05
            df_statstable['sig_color'] = 'k'
            df_statstable['sig_label'] = ' '
            # df_statstable['sig_color'] = df_statstable['p_value'].apply(
            #                                     lambda x: 'r' if x > use_pvalue else 'k')
            # # df_statstable['sig_label'] = df_statstable['p_value'].apply(
            #                                     lambda x: 'not significant' if x > use_pvalue else 'significant')

            # Create the map only once for both lines ----------------------------------------------
            if iscenario == 1:
                # get start middle and ends of lines for labelling
                start_end_pts = gdf_map[gdf_map[transect_column] == 'Strip'].iloc[[0, -1]]
                middle_pts = gdf_map[gdf_map['PointID'] == gdf_map.iloc[int(len(gdf_map) / 2)]['PointID']]

                # create lines

                # Aggregate points with the GroupBy - for shapely 1.3.2+ use the point objects
                try:
                    gdf_lines = gdf_map.groupby(['Strip_Name'])['geometry'].apply(
                        lambda x: LineString(x.tolist()) if x.size > 1 else x.tolist())
                except ValueError:
                    gdf_lines = gdf_map.groupby(['Strip_Name'])['geometry'].apply(
                        lambda x: LineString([(p.x, p.y) for p in x]))

                gdf_lines = GeoDataFrame(gdf_lines, geometry='geometry')

                fig_map, ax = plt.subplots(figsize=(5, 5))

                # fix the extent to the vector
                ax.set_aspect('equal')
                ax.set_title('Map for Trial {}'.format(line_id), fontsize=8)
                ax.set_xlim(left=gdf_strip.total_bounds[0] - 50, right=gdf_strip.total_bounds[2] + 50)
                ax.set_ylim(bottom=gdf_strip.total_bounds[1] - 50, top=gdf_strip.total_bounds[3] + 50)

                # plot the line
                gdf_lines.plot(ax=ax, cmap='rainbow', legend=True)
                # _ = gdf_map.plot(ax=ax,marker='o',column='Strip_Name',cmap='rainbow',markersize=5)

                ax.tick_params(which='major', width=0.75, length=2.5, labelsize=4)
                # ,top=True,labeltop=True,bottom=False,labelbottom=False)

                plt.setp(ax.spines.values(), linewidth=0.5)

                ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{:,.0f}".format(x)))
                ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{:,.0f}".format(x)))
                # plt.xticks(rotation='vertical')

                anno_kwargs = {'size'        : 7, 'va': 'center', 'ha': 'center',
                               'path_effects': [pe.withStroke(linewidth=4, foreground="white")]}

                middle_pts.apply(lambda x: ax.annotate(x['Strip_Name'],
                                                       xy=(x.geometry.x, x.geometry.y),
                                                       rotation=x['label_angle'], **anno_kwargs),
                                 axis=1)

                ax.annotate('Start', xy=(start_end_pts.iloc[0].geometry.x,
                                         start_end_pts.iloc[0].geometry.y),
                            rotation=start_end_pts.iloc[0]['label_angle'], **anno_kwargs)

                ax.annotate('End', xy=(start_end_pts.iloc[-1].geometry.x,
                                       start_end_pts.iloc[-1].geometry.y),
                            rotation=start_end_pts.iloc[-1]['label_angle'], **anno_kwargs)

                map_file = os.path.join(out_folder, "Trial-{}_map.png".format(line_id))

                plt.savefig(map_file, dpi=300)
                plt.close()

                if config.get_debug_mode():
                    LOGGER.info('{:<30}\t{:>10}   {dur:<15} {}'.format(
                        'Map Saved', '', map_file, dur=str(timedelta(seconds=time.time() - loop_time))))

            # Create the graphs ----------------------------------------------
            fig_graph = plt.figure(figsize=(15, 10))  # width, height
            fig_graph.subplots_adjust(hspace=0, wspace=0.05)
            gs = gridspec.GridSpec(3, 1, height_ratios=[3, 3, 1])

            if config.get_debug_mode():
                gs = gridspec.GridSpec(4, 1, height_ratios=[3, 3, 1, 3])

            axs = []
            for igs, ss in enumerate(gs):
                if igs == 0:
                    axs.append(fig_graph.add_subplot(ss))
                else:
                    axs.append(fig_graph.add_subplot(ss, sharex=axs[0]))

            # add a title
            axs[0].set_title('Strip trial analysis for Trial {} - {} strip'.format(line_id, offset_subst), fontsize=16)

            df_statstable.plot(x='DistOnLine', y=column_names['Value'][0], marker='.', ax=axs[0], label='Treatment')

            df_statstable.plot(x='DistOnLine', y=control_mean, marker='.', ax=axs[0], label='Control')

            for name, group in df_statstable.groupby(['Strip Zone', 'zone_marker', 'sig_color', 'sig_label']):
                group.plot.scatter(x='DistOnLine', y='av_treat_dif', marker=name[1], s=25,
                                   c=name[2], ax=axs[1], zorder=2,
                                   label='Zone {} {}'.format(name[0], name[3]))  # ,s=25

                group.plot.scatter(x='DistOnLine', y='RI', marker=name[1], s=25, c='k',
                                   ax=axs[2], zorder=2,
                                   label='Zone {}'.format(name[0]))  # s=25

                # if config.get_debug_mode():
                #     group.plot.scatter(x='DistOnLine', y='p_value', marker=name[1], s=25, c='k',
                #                        ax=axs[3], label='Zone {}'.format(name[0]), zorder=2)

            axs[0].set(ylabel="Treatment Units")
            axs[1].set(ylabel="Treatment Difference")
            axs[2].set(ylabel="Response\nIndex")

            # if config.get_debug_mode():
            #     axs[3].set(ylabel="p value")
            #     for i, ea in enumerate([0.05, 0.01, 0.005, 0.001]):
            #         axs[3].axhline(y=ea, color=colour_list[i], linestyle='-', alpha=0.8, zorder=1,
            #                        label='p={}'.format(ea))

            for iax, ea_ax in enumerate(axs):
                ea_ax.grid(True, which='major', axis='x')
                ea_ax.set(xlabel="Distance (m) from start of strip (see map)")
                # increase graph outline
                plt.setp(ea_ax.spines.values(), linewidth=1)

                # Shrink current axis by 20%  https://stackoverflow.com/a/4701285
                box = ea_ax.get_position()
                ea_ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ea_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

                # remove duplicates labels from legend  https://stackoverflow.com/a/13589144
                handles, labels = ea_ax.get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))

                if config.get_debug_mode() and iax == 1:
                    ea_ax.legend(by_label.values(), by_label.keys(), loc='center left',
                                 bbox_to_anchor=(1, 0.5), edgecolor='w')
                    # title='Using a p-value of {} '.format(use_pvalue))
                else:
                    ea_ax.legend(by_label.values(), by_label.keys(), loc='center left',
                                 bbox_to_anchor=(1, 0.5), edgecolor='w')

            plt.savefig(file_path_noext + '_graph.png')
            plt.close()
            if config.get_debug_mode():
                LOGGER.info('{:<30}\t{:>10}   {dur:<15} {}'.format(
                    'Saved graph', '', file_path_noext + '_graph.png',
                    dur=str(timedelta(seconds=time.time() - loop_time))))

        LOGGER.info('{:<30} {:>15} {dur}'.format(
            'T-Test for Line {} Completed'.format(line_id), status,
            dur=str(timedelta(seconds=time.time() - loop_time))))

    return df_table


def persistor_all_years(raster_files, output_tif, greater_than, target_percentage):
    """Determine the performance persistence of yield by across multiple years as described in
    Bramley and Hamilton (2005)

    The "Target over all years" method assigns a value to each pixel to indicate the number of
    instances (in the raster list) in which that pixel was either less than or greater than
    the mean (+/- a nominated percentage) of that raster.

     All input rasters MUST overlap and have the same coordinate system and pixel size.

     If a path is omitted from output_tif it will be created in your temp folder.

    References:
        Bramley RGV, Hamilton RP (2005) Understanding variability in winegrape production systems
        1. Within vineyard variation in yield over several vintages. Australian Journal Of Grape And
        Wine Research 10, 32-45. doi:10.1111/j.1755-0238.2004.tb00006.x.

    Args:
        raster_files (List[str]): List of rasters to use as inputs
        output_tif (str): Output TIF file
        greater_than (bool): if true test above (gt) the mean
        target_percentage (int): the percent variation either above/below the mean.
                   This should be a integer between -50 and 50

    Returns:
        str: The output tif name
    """

    if not isinstance(target_percentage, int):
        raise TypeError('target_percentage must an integer between -50 to 50.')

    if target_percentage < -50 or target_percentage > 50:
        raise ValueError('target_percentage must an integer between -50 to 50.')

    if not isinstance(greater_than, (int, bool)):
        raise TypeError('greater_than must be boolean.')

    if output_tif is not None and not os.path.isabs(output_tif):
        output_tif = os.path.join(TEMPDIR, output_tif)

    if not isinstance(raster_files, list):
        raise TypeError('Invalid Type: raster_files should be a list')

    if len(raster_files) == 0:
        raise TypeError('Invalid Type: Empty list of raster files')

    start_time = time.time()

    # combine and align grid and slip to common area
    image_combined, _ = stack_and_clip_rasters(raster_files, use_common=True,
                                               output_tif='allyrs_combined.tif')

    with rasterio.open(image_combined, 'r') as src:
        kwargs = src.meta.copy()

        for i, band in enumerate(src.read(masked=True), start=1):
            cutoff = np.nanmean(band) + (np.nanmean(band) * (target_percentage / 100.00))
            band_new = band.copy()

            # classify bands
            if greater_than:
                band_new = np.where(band > cutoff, 1, 0)
            else:
                band_new = np.where(band < cutoff, 1, 0)

            # reapply nodata values
            band_new[np.where(band.mask)] = -9999
            band_new = band_new.astype(np.int16)

            # cumulative sum of bands
            if i == 1:
                sum_of_rasters = band_new
            else:
                sum_of_rasters = sum_of_rasters + band_new

            sum_of_rasters[band.mask] = -9999

            del band, band_new

    kwargs.pop('count')
    kwargs.update({'nodata': -9999, 'dtype': np.int16})
    with rasterio.open(output_tif, 'w', count=1, **kwargs) as dst_allyrs:
        dst_allyrs.write(sum_of_rasters.astype(np.int16), 1)

    arg_str = '{} {}%'.format('>' if greater_than else '<', target_percentage)
    LOGGER.info('{:<30} {:>15} {dur}'.format('Persistor All Years Completed', arg_str,
                                             dur=str(timedelta(seconds=time.time() - start_time))))

    # cleanup Temp files
    if not config.get_debug_mode() and TEMPDIR in image_combined:
        os.remove(image_combined)

    return output_tif


def persistor_target_probability(upper_raster_files, upper_percentage, upper_probability,
                                 lower_raster_files, lower_percentage, lower_probability,
                                 output_tif):
    """Determine the probability of a performance being exceeded or not being met with an upper and
     lower limit as described in Bramley and Hamilton (2005).

    The "Target probability" method builds on the target over all years method, in that it includes
    an upper range (i.e. cells with a given frequency of values that are above the mean +/- a
     given percentage) and a lower range (i.e. cells with a given frequency of values that are
     below the mean +/- a given percentage).

    A value is assigned to each pixel which indicates whether the performance in that pixel over
    a given proportion of years is:
        a)	Greater than the mean plus or minus the nominated percentage (value = 1)
        b)	Less than the mean plus or minus the nominated percentage (value = -1)
        The remaining pixels which do not fall into category a) or b) are given a value of 0.

    All input rasters MUST overlap and have the same coordinate system and pixel size.

    If a path is omitted from output_tif it will be created in your temp folder.

    References:
        Bramley RGV, Hamilton RP (2005) Understanding variability in winegrape production systems
        1. Within vineyard variation in yield over several vintages. Australian Journal Of Grape
        And Wine Research 10, 32-45. doi:10.1111/j.1755-0238.2004.tb00006.x.

    Args:
        upper_raster_files (List[str]): List of rasters to used for the analysis of the
                    UPPER category
        upper_percentage (int): the percent variation either above/below the mean to apply
                    to the UPPER raster category.
        upper_probability (int):the probability percentage to apply to the LOWER category

        lower_raster_files (List[str]): List of rasters to used for the analysis of the
                    LOWER category
        lower_percentage (int): the percent variation either above/below the mean to apply
                    to the LOWER raster category.
        lower_probability (int): the probability percentage to apply to the LOWER category

        output_tif (str): Output TIF file

    """

    for arg_check in [('upper_raster_files', upper_raster_files),
                      ('lower_raster_files', lower_raster_files)]:
        if not isinstance(arg_check[-1], list):
            raise TypeError('Invalid Type: {} should be a list'.format(arg_check[0]))

        if len(arg_check[-1]) == 0:
            raise TypeError('Invalid Type: {} is a empty list'.format(arg_check[0]))

    for arg_check in [('upper_percentage', upper_percentage),
                      ('upper_probability', upper_probability),
                      ('lower_percentage', lower_percentage),
                      ('lower_probability', lower_probability)]:
        if not isinstance(arg_check[1], int):
            raise TypeError('{} must an integer'.format(arg_check[0]))

    for arg_check in [('upper_percentage', upper_percentage),
                      ('lower_percentage', lower_percentage)]:
        if not isinstance(arg_check[1], int):
            raise ValueError('{} must an integer between -50 to 50.'.format(arg_check[0]))

    if output_tif is None or output_tif == '':
        raise TypeError('Please specify an output filename')

    if not os.path.exists(os.path.dirname(output_tif)):
        raise IOError('Output directory {} does not exist'.format(os.path.dirname(output_tif)))

    start_time = time.time()
    out_prefix, out_ext = os.path.splitext(os.path.basename(output_tif))

    # Process upper ------------------------------------------------------------------------------
    temp_file_list = [os.path.join(TEMPDIR, out_prefix + '_1_upperyears' + out_ext)]
    upper_tif = persistor_all_years(upper_raster_files,
                                    temp_file_list[-1],
                                    greater_than=True,
                                    target_percentage=upper_percentage)

    upper_cutoff = (upper_probability / 100.00) * len(upper_raster_files)
    with rasterio.open(upper_tif, 'r') as src:
        kwargs = src.meta.copy()
        band = src.read(masked=True)

        # Greater than test returns boolean which can be converted to 0,1
        data_u = np.where(~band.mask, (band >= upper_cutoff).astype(np.int16), -9999)

    if config.get_debug_mode():
        temp_file_list += [os.path.join(TEMPDIR, out_prefix + '_2gt_upperprob' + out_ext)]
        with rasterio.open(temp_file_list[-1], 'w', **kwargs) as dst:
            dst.write(data_u)

    del band

    # Process lower ------------------------------------------------------------------------------
    temp_file_list += [os.path.join(TEMPDIR, out_prefix + '_1_loweryears' + out_ext)]
    lower_tif = persistor_all_years(lower_raster_files,
                                    temp_file_list[-1],
                                    greater_than=False,
                                    target_percentage=lower_percentage)

    lower_cutoff = (lower_probability / 100.00) * len(lower_raster_files)
    with rasterio.open(lower_tif, 'r') as src:
        band = src.read(masked=True)

        # (band >= lower_cutoff) returns booleans (0,1) we need (0,-1) so
        # convert to integer and use np.negative to switch it while maintaining the nodata values.
        data_l = np.where(~band.mask, np.negative((band >= lower_cutoff).astype(np.int16)), -9999)

        if config.get_debug_mode():
            temp_file_list += [os.path.join(TEMPDIR, out_prefix + '_2gtn_lowerprob' + out_ext)]
            with rasterio.open(temp_file_list[-1], 'w', **kwargs) as dst:
                dst.write(data_l)

        with rasterio.open(output_tif, 'w', **kwargs) as dst:
            result = np.where(~band.mask, data_l + data_u, -9999)
            dst.write(result.astype(np.int16))

    LOGGER.info('{:<30} {:>15} {dur}'.format('Persistor Target Probability Completed', '',
                                             dur=str(timedelta(seconds=time.time() - start_time))))

    # clean up of intermediate files
    if len(temp_file_list) > 0 and not config.get_debug_mode():
        for ea in temp_file_list:
            for f in glob.glob(os.path.splitext(ea)[0] + '*'):
                os.remove(f)
    return
