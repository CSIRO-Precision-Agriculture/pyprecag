import logging
import os
import six
import warnings

from packaging.version import parse as parse_version

from osgeo import osr, gdal
import pyproj

from pyproj import CRS
import rasterio as rio
import geopandas as gpd

from ._compat import SHAPELY_GE_20

if SHAPELY_GE_20:
    from shapely import Point, box
else:
    from shapely.geometry import Point, box

from . import config

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())

class crs:
    def __init__(self):
        warnings.warn("pyprecag.crs has been deprecated in favor of `pyproj`, `geopandas.crs` or `rasterio.crs` "
                      "and will be removed in a future release", PendingDeprecationWarning, 2)

        self.srs = None  # Spatial Reference Object
        self.crs_wkt = None
        self.proj4 = None
        self.epsg_number = None
        self.epsg = None
        self.epsg_predicted = False  # Did the epsg_number get set via the online lookup.
        self.pyproj_version = pyproj.__version__

    def set_epsg(self, code):
        """Given an integer code, set the epsg_number and the EPSG-like mapping.

         Note: the input code is not validated against an EPSG database.
         Args:
            code (int): An integer representing the EPSG code

         """
        if code is None or code == 0:
            warnings.warn('EPSG Code is None or 0 and cant be set.')
            return

        if isinstance(code, six.string_types):
            try:
                code = int(code)
            except:
                raise TypeError('EPSG code ({}) must be a positive integer'.format(code))

        if not isinstance(code, int):
            raise TypeError('EPSG Code ({}) must be a positive integer'.format(code))
        if code <= 0:
            raise ValueError('EPSG Code ({}) must be a positive integer above 0'.format(code))

        self.epsg = from_epsg(code)
        self.epsg_number = code

    def getFromWKT(self, crs_wkt, bUpdateByEPSG=False):
        """ Get crs attributes via a coordinate systems WKT.
            Used for gathering coordinate system attributes from fiona.crs.crs_wkt.
            If the EPSG number can be determined, bUpdateByEPSG provides the option to update the proj4, crs_wkt and
            SRS attributes to match the official EPSG database attributes
        Args:
            crs_wkt (unicode):    A string representing the coordinate system well known text
            bUpdateByEPSG (bool): If True and the epsg_number is predicted, overwrite crs attributes with those from
                                  official epsg_number database.

        """
        if crs_wkt == '':
            warnings.warn('crs_wkt is blank. No coordinate information set')
            return

        source = osr.SpatialReference()
        source.ImportFromWkt(crs_wkt)

        if source is not None:
            self.crs_wkt = crs_wkt
            self.srs = source
            self.proj4 = source.ExportToProj4()

            source.AutoIdentifyEPSG()
            code = source.GetAuthorityCode(None)

            if code is None:
                if parse_version(pyproj.__version__) >= parse_version('2.4.0'):
                    pyproj_crs = pyproj.CRS(source.GetName())
                    code = pyproj_crs.to_epsg()
            if code is None:
                code = self.getEPSGFromSRS(source, bUpdateToCorrectDefn=bUpdateByEPSG)

            if code:
                self.set_epsg(int(code))

    def getFromEPSG(self, epsg):
        """Create OGR Spatial Reference Object for an epsg_number number
        Args:
            epsg (int):    A Valid EPSG Number

        Returns:
            osgeo.osr.SpatialReference: The Spatial Reference Object for the input EPSG

        """
        if isinstance(epsg, six.string_types):
            try:
                epsg = int(epsg.replace('EPSG:', ''))
            except:
                raise TypeError('EPSG must be a number. Got - {}'.format(epsg))

        if not isinstance(epsg, six.integer_types):
            raise TypeError('EPSG must be a number. Got - {}'.format(epsg))

        if epsg is None or epsg == 0:
            return None

        srs_obj = osr.SpatialReference()
        srs_obj.ImportFromEPSG(epsg)

        if srs_obj is not None:
            self.srs = srs_obj
            self.crs_wkt = srs_obj.ExportToWkt()
            self.proj4 = srs_obj.ExportToProj4()
            self.set_epsg(epsg)
            self.epsg_predicted = False

    def getEPSGFromSRS(self, osr_srs, bUpdateToCorrectDefn=False):
        """Get the EPSG number for a Spatial Reference System.

           If the Spatial reference system does not contain an EPSG number it will attempt to find one.

           adapted from : https://gis.stackexchange.com/a/8888

        Args:
            osr_srs (osgeo.osr.SpatialReference):
            bUpdateToCorrectDefn (bool): overwrite existing with the found definition
        Returns:
            int: EPSG Number for the matched Spatial Reference System.

        """

        if osr_srs is None:
            return
        orig_srs = osr_srs.Clone()

        for i in [0, 1]:
            osr_srs.AutoIdentifyEPSG()
            epsg = osr_srs.GetAuthorityCode(None)
            if epsg is not None:
                try:
                    epsg = int(epsg)
                except:
                    epsg = None

                if epsg > 0:
                    return epsg

            # For the second loop try the after using MorphFromESRI
            osr_srs.MorphFromESRI()

        # reset back to original to undo MorphFromESRI
        osr_srs = orig_srs.Clone()

        # Then through a lookup
        if epsg is None or epsg == 0:
            self.epsg_predicted = True
            if osr_srs.IsProjected():
                srsName = osr_srs.GetAttrValue("PROJCS", 0)
            else:
                srsName = osr_srs.GetAttrValue("GEOGCS", 0)

            self.epsg_predicted = True

            # this is copied from QGIS 3.10.3 as a temporary solution
            crs_database = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'proj.db')

            import sqlite3
            import pandas as pd

            connection = sqlite3.connect(crs_database)
            cursor = connection.cursor()
            cs_df = pd.read_sql_query('Select  * from  alias_name', connection)
            del cursor
            connection.close()

            cs_df = cs_df[cs_df['alt_name'].str.contains(srsName.replace('GDA', ''), case=False)]

            if len(cs_df) == 1:
                epsg = int(cs_df['code'])

            if len(cs_df):
                LOGGER.debug('Matched to ESRI alias {}'.format(srsName))

        if bUpdateToCorrectDefn and epsg > 0:
            self.getFromEPSG(int(epsg))
            self.set_epsg(epsg)

            return int(epsg)


def from_epsg(epsg_number):
    """
    Get a coordinate system for the epsg number.
    if pyproj is pre 2.4.0 then use a string created using fiona.crs.from_epsg, otherwise use the pyproj option.
    Args:
        epsg_number: The epsg number representing the coordinate system.

    Returns:
        object: pyproj.crs     The newer pyproj 2.4+ crs object compatible with Proj 6+
                or
               string          Derived from epsg using fiona.
    """
    warnings.warn('pyprecag.crs.from_epsg() in favor of `pyproj, geopandas or rasterio solutions` '
                  'see https://pyproj4.github.io/pyproj/stable/crs_compatibility.html',
                  PendingDeprecationWarning, stacklevel=2)

    if parse_version(pyproj.__version__) >= parse_version('2.4.0'):
        crs = pyproj.CRS.from_epsg(epsg_number)
    else:
        from fiona import crs as fio_crs
        crs = fio_crs.from_epsg(epsg_number)
    return crs


def getCRSfromRasterFile(raster_file):
    """Create a CRS object from a raster file.

    Args:
        raster_file (str): The path and filename to a raster file

    Returns:
        pyprecag.crs.crs: An object representing the raster file's coordinate system.

    """
    warnings.warn('pyprecag.crs.getCRSfromRasterFile() is deprecated in favor of `rasterio`'
                  ' - see https://rasterio.readthedocs.io/en/stable/api/rasterio.crs.html',
                  PendingDeprecationWarning, stacklevel=2)

    gdalRaster = gdal.Open(os.path.normpath(raster_file))
    rast_crs = crs()
    rast_crs.getFromWKT(gdalRaster.GetProjectionRef())
    del gdalRaster
    return rast_crs


def getUTMfromWGS84(longitude, latitude):
    """ Calculate the UTM Zonal projection from a set of lats&longs

    This is useful for converting  decimal degrees to metres for area perimeter etc.
    utm for WGS84 is selected as it is global, and wgs84 is commonly used as the GPS coordinate system.

    It will return the zone, and the spatial reference objects for the utm & wgs84 projections required
    for coordinate or feature transformation.

    Args:
        longitude (float): a floating number representing longitude
        latitude (float): a floating number representing latitude

    Returns:
        Tuple[int, osgeo.osr.SpatialReference, osgeo.osr.SpatialReference]: UTM Zone, utm SRS, WGS84 SRS

    """

    # Based On :https://stackoverflow.com/a/10239676
    def get_utm_zone(longitude):
        return int(1 + (longitude + 180.0) / 6.0)

    zone = get_utm_zone(longitude)

    # Determines if given latitude is a northern for UTM
    is_northern = int(latitude > 0.0)

    osr_srs = osr.SpatialReference()
    osr_srs.SetWellKnownGeogCS("WGS84")  # Set geographic coordinate system to handle latitude/longitude
    osr_srs.SetUTM(zone, is_northern)
    osr_srs.AutoIdentifyEPSG()
    out_epsg = int(osr_srs.GetAuthorityCode(None))
    osr_srs.ImportFromEPSG(out_epsg)
    wgs84_crs = osr_srs.CloneGeogCS()  # Clone ONLY the geographic coordinate system

    return zone, osr_srs, wgs84_crs, out_epsg


def match_crs_to_epsg(input_crs):
    """
    when an epsg for a coordinate system can't be determined loop through different
    WKT versions and try and find a match.
    Args:
        input_crs (pyproj/rasterio/geopandas CRS): the input coordinate system to search on.
                                                   this could be from pyproj, rasterio or pandas.

    Returns:
        pyproj.crs.crs.CRS:
    """

    if not isinstance(input_crs, (pyproj.crs.CRS, int, rio.crs.CRS)):
        raise TypeError(f'input_crs must a pyproj CRS object or epsg number. Got - {type(input_crs)}')

    if isinstance(input_crs, rio.crs.CRS):
        input_crs = CRS.from_user_input(input_crs)

    for wkt_version in pyproj.enums.WktVersion:
        wkt = input_crs.to_wkt(wkt_version)
        input_crs = CRS.from_user_input(wkt)
        if input_crs.to_epsg():
            return CRS.from_epsg(input_crs.to_epsg())

    warnings.warn('Could not match possibly due to a custom or poorly defined coordinate system')
    return None


def get_projected_epsg_for_bounds(bounds, bounds_crs):
    """"
        Determine a projected coordinate system from a bounding box (xmin, ymin, xmax,ymax). This allows for local
        projected coordinate system to be preferred over the standard UTM zones. The local coordinate system are
        listed in the pyprecag/config.json configuration file.

    Args:
        bounds(list or tuple): [xmin, ymin, xmax,ymax]
        bounds_crs (pyproj/raster CRS object or epsg):  a pyproj/raster coordinate system or epsg for the bounding box.

    Returns:
        int: projected coordinate system as an EPSG number.

    """
    if not bounds_crs or bounds_crs == 0:
        return

    if not isinstance(bounds_crs, (pyproj.crs.CRS, int, rio.crs.CRS)):
        raise TypeError('Needs a pyproj/rasterio/geopandas CRS object or epsg number')
        return

    crs_type = type(bounds_crs)

    if isinstance(bounds_crs, rio.crs.CRS):
        bounds_crs = CRS.from_user_input(bounds_crs)
    elif isinstance(bounds_crs, int):
        bounds_crs = CRS.from_epsg(bounds_crs)

    # check for a properly formed crs
    if not bounds_crs.to_epsg():
        bounds_crs = match_crs_to_epsg(bounds_crs)

    if bounds_crs.is_projected:
        return bounds_crs.to_epsg()

    epsg = get_projected_epsg_for_point(box(*bounds).representative_point(), bounds_crs)

    return epsg


def get_projected_epsg_for_point(point, point_crs):
    """
        Determine a projected coordinate system for a shapely point. This allows for local projected coordinate system
        to be preferred over the standard UTM zones. The local coordinate system are listed
        in the pyprecag/config.json configuration file.
    Args:
        point (shapely.geometry.point.Point):
        point_crs (pyproj.crs.crs.CRS):

    Returns:
        int:  EPSG for the projected coordinate system
    """
    if not point_crs or point_crs == 0:
        return

    if not isinstance(point_crs, (pyproj.crs.CRS, int, rio.crs.CRS)):
        raise TypeError('Needs a pyproj/rasterio/geopandas CRS object or epsg number')
        return

    crs_type = type(point_crs)

    if isinstance(point_crs, rio.crs.CRS):
        point_crs = CRS.from_user_input(point_crs)
    elif isinstance(point_crs, int):
        point_crs = CRS.from_epsg(point_crs)

    # check for a properly formed crs
    if not point_crs.to_epsg():
        point_crs = match_crs_to_epsg(point_crs)

    if point_crs.is_projected:
        return point_crs.to_epsg()

    # create bounds for local coordinate systems ie not global UTM
    rec = []
    for epsg in config.get_config_key('local_projected_epsg'):
        l_crs = CRS.from_epsg(epsg)
        rec.append({'Name': l_crs.name, 'Datum': l_crs.datum.name, 'EPSG': epsg,
                    'geometry': box(*l_crs.area_of_use.bounds)})

    crs_gdf = gpd.GeoDataFrame(rec, geometry='geometry', crs=4326)

    # create point from xy and reproject to 4326
    pt_gdf = gpd.GeoDataFrame(geometry=[point], crs=point_crs).to_crs(crs_gdf.crs)

    # find out if our target is within one local crs extents and if not use global utm zones
    result = gpd.sjoin(left_df=crs_gdf, right_df=pt_gdf.to_crs(crs_gdf.crs), how='inner', predicate='intersects')

    if len(result) == 1:
        epsg = result.at[result.index[0], 'EPSG']
    else:
        epsg = pt_gdf.estimate_utm_crs().to_epsg()

    return epsg


def getProjectedCRSForXY(x_coord, y_coord, xy_epsg=4326):
    """ Calculate the Zonal projected coordinate system from a set of xy coordinates.

        Coordinates will be reprojected to geographics wgs84(4326) if required to identify the correct
        projected coordinate system.

        utm for WGS84 is selected as it is global, and wgs84 is commonly used as the GPS coordinate system.

        If input coordinates are within Australia the result will be in Map Grid of Australia(MGA) GDA 1994 coordinate
        system, otherwise the appropriate zonal utm WGS84 projected coordinate system will be used.

        This is useful for
            - converting decimal degrees to metres for area perimeter etc.
            - converting sUTM to GDA mga

        Args:
            x_coord (float): A floating number representing an easting, x or longitude
            y_coord (float): A floating number representing a northing, y or latitude
            xy_epsg (int): The epsg_number for the x_coord & y_coord coordinates.
                               This could geographic(WGS84, or GDA94) or projected(UTM or MGA GDA)

        Returns:
            pyprecag.crs : A pyprecag.crs object defining the WGS84 Zone Projected Coordinate System

        """

    # Coordinates need to be in wgs84 so project them
    if xy_epsg != 4326:
        if parse_version(pyproj.__version__) >= parse_version('2.4.0'):
            inProj = from_epsg(xy_epsg)  # Proj(init='epsg:{}'.format(xy_epsg))
            outProj = from_epsg(4326)  # Proj(init='epsg:4326')

            from pyproj import Transformer
            transformer = Transformer.from_crs(inProj, outProj, always_xy=True)
            longitude, latitude = transformer.transform(x_coord, y_coord)
        else:
            # Based On :https://stackoverflow.com/a/10239676
            from pyproj import Proj, transform
            inProj = Proj(init='epsg:{}'.format(xy_epsg))
            outProj = Proj(init='epsg:4326')
            longitude, latitude = transform(inProj, outProj, x_coord, y_coord)

    else:
        longitude, latitude = x_coord, y_coord

    utm_zone, _, _, out_epsg = getUTMfromWGS84(longitude, latitude)

    osr_srs = osr.SpatialReference()
    # Use Australian or New Zealand local systems otherwise use global wgs84 UTM zones.
    if (108.0 <= longitude <= 155.0) and (-45.0 <= latitude <= -10.0):
        # set to Australian GDA94 MGA Zone XX
        out_epsg = int('283{}'.format(utm_zone))
    elif (166.33 <= longitude <= 178.6) and (-47.4 <= latitude <= -34.0):
        # set to NZGD2000 / New Zealand Transverse Mercator 2000
        out_epsg = 2193

    utmzone_crs = crs()
    utmzone_crs.getFromEPSG(out_epsg)

    if utmzone_crs.srs.IsProjected():
        return utmzone_crs
    else:
        return None


def distance_metres_to_dd(longitude, latitude, distance_metres):
    """ Converts a distance in metres to decimal degrees. It presumes the longitude/lats are in WGS84.

    Workflow -> using a Long/Lat, reproject to UTM
             -> Add distance to the easting
             -> reproject back to WGS84.
             -> Calculate distance between the two sets of coordinates.

    Args:
        longitude (float): a floating number representing longitude
        latitude (float): a floating number representing latitude
        distance_metres (float): a distance in metres to convert.

    Returns:
        float: the distance in decimal degrees.

    """

    for argCheck in [('longitude', longitude), ('latitude', latitude)]:
        if not isinstance(argCheck[1], float):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    if not isinstance(distance_metres, six.integer_types + (float,)):
        raise TypeError('distance_metres must be a floating number.')

    # get the required Spatial reference systems
    utm_crs = getProjectedCRSForXY(longitude, latitude)
    wgs84SRS = utm_crs.srs.CloneGeogCS()
    # zone, utmSRS, wgs84SRS = getUTMfromWGS84(longitude, latitude)

    # create transform component
    wgs842utm_transform = osr.CoordinateTransformation(wgs84SRS, utm_crs.srs)  # (<from>, <to>)

    # create transform component
    utm2wgs84_transform = osr.CoordinateTransformation(utm_crs.srs, wgs84SRS)  # (<from>, <to>)

    # Project to UTM
    easting, northing, alt = wgs842utm_transform.TransformPoint(longitude, latitude)

    # if easting/northing return infity, then swap lat and long over.
    if easting == float("inf"):
        easting, northing, alt = wgs842utm_transform.TransformPoint(latitude, longitude, 0)

        # Project back to WGS84 after adding distance to eastings. returns easting, northing, altitude
        newLat, newLong, alt = utm2wgs84_transform.TransformPoint(easting + distance_metres, northing, 0)
    else:
        # Project back to WGS84 after adding distance to eastings. returns easting, northing, altitude
        newLong, newLat, alt = utm2wgs84_transform.TransformPoint(easting + distance_metres, northing, 0)

    # Create a point objects.
    origPoint = Point(longitude, latitude)
    adjPoint = Point(newLong, newLat)

    # Calculate Distance
    distance_dd = origPoint.distance(adjPoint)

    # if the input was negative make sure the output is too
    if distance_metres < 0:
        distance_dd = -distance_dd

    return distance_dd
