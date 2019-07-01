Changelog
=========
**0.3.0 (2019-06-30)**
 New Tools
  * #24 Persistor.  
  * #22 *t*-test analysis of strip trials.  
 Enhancements
  * #30 Low spatial density kriging using VESPER.  
  * #22 Updated processing.create_points_along_line to allow for line offset to use compass points instead of Offset 1 etc.
  
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
