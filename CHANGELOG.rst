Changelog
=========
** 0.4.2 (2023-10-01) **
 + Now supports rasterio 1.3.8, geopandas 0.13.2, GDAL 3.7.1
 + Improved method for creating strip trial vectors
 + Faster thin by distance method.
 + implemented pyproj6/geopandas/rasterio crs functionality and deprecated old functions.
 Deprecated functions
  + ``convert.drop_z()`` replaced by ``shapely.force_2d()``
  + ``vector_ops.explode_multi_part_features()`` replaced by ``GeoDataFrame.explode()``
  + ``vector_ops.calculate_area_length_in_metres()`` replaced by ``GeoDataFrame.to_crs(GeoDataFrame.estimate_utm_crs().to_epsg()).geometry.area`` or equivalent
 Bug Fixes
  + #49 writing vesper grid is slow
  + #50 fixes for Dataframe.append() deprecation in pandas/geopandas.

** 0.4.1 (2023-02-07)**
 + Support for newer versions of rasterio 1.3.3, geopandas 0.11.1, GDAL 3.6.1
 + #45 Fix to k-means clustering
 + #44 Fix to clean/trim
 + #41 Fix for empty shapefile in Clean/trim
 + Added vesper statistics as tags in output prediction tif

** 0.4.0 (2020-08-07)**
 * #40 Port to Python3
 * Updated Blockgrid to support a group-by batch mode.
 * ``block_size`` argument removed from ``kriging_ops.prepare_for_vesper_krige``

**0.3.1 (2020-05-27)**
 Bug Fixes
  * Support for Geopandas (upto 0.5.1)
  * Fix coordinate system lookup error (PAT Issue `#42 <https://github.com/CSIRO-Precision-Agriculture/PAT_QGIS_Plugin/issues/42>`_ )

**0.3.0 (2019-06-30)**
 New Tools
  * #24 Persistor.
  * #22 *t*-test analysis of strip trials.
 Enhancements
  * #30 Low spatial density kriging using VESPER.  
  * #22 Updated ``processing.create_points_along_line`` to allow for line offset to use compass points instead of Offset 1 etc.
 Minor Bug Fixes.

**0.2.2 (2019-02-27)**
  #15 Changes to processing.create_points_along_line.
   * Now works with geopandas 0.3.0 and 0.4.0.
   * Renamed columns names: side -> Transect, line_dist -> DistOnLine.
   * Change Transect/side attributes: C -> Centre, L -> Offset 1, R -> Offset 2.

**0.2.1 (2019-01-30)**  
 *  #15 Create points along lines now correctly saves the lines to shapefile.

**0.2.0 (2019-01-24)**
 * New Feature - Create points for a strip trial along a central line and offset left and right at distance.

**0.1.1 (2019-01-22)**
 * New Feature - Create zones using *k*-means clustering.

**0.0.4 (2018-12-20)**
 * Initial load into GitHub.
