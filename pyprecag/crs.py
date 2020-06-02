import json
import logging

import os
import socket
import warnings
from urllib import urlencode
import urllib2

from fiona.crs import from_string, from_epsg
from osgeo import osr, gdal
from shapely import geometry

from . import config
from .errors import SpatialReferenceError

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
# LOGGER.setLevel("DEBUG")
# DEBUG = config.get_debug_mode()  # LOGGER.isEnabledFor(logging.DEBUG)


class crs:
    def __init__(self):
        self.srs = None  # Spatial Reference Object
        self.crs_wkt = None
        self.proj4 = None
        self.epsg_number = None
        self.epsg = None
        self.epsg_predicted = False  # Did the epsg_number get set via the online lookup.

    def set_epsg(self, code):
        """Given an integer code, set the epsg_number and the EPSG-like mapping.

         Note: the input code is not validated against an EPSG database.
         Args:
            code (int): An integer representing the EPSG code

         """
        if code is None:
            warnings.warn('EPSG Code is None and cant be set.')
            return

        if isinstance(code, str):
            try:
                code = int(code)
            except:
                raise TypeError('EPSG code must be a positive integer')

        if not isinstance(code, int):
            raise TypeError('EPSG Code must be a positive integer')
        if code <= 0:
            raise ValueError('EPSG Code must be a positive integer')

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
                code = self.getEPSGFromSRS(source,  bUpdateToCorrectDefn=bUpdateByEPSG)
            else:
                code = int(source.GetAuthorityCode(None))

            if code is not None:
                self.set_epsg(code)

    def getFromEPSG(self, epsg):
        """Create OGR Spatial Reference Object for an epsg_number number
        Args:
            epsg (int):    A Valid EPSG Number

        Returns:
            osgeo.osr.SpatialReference: The Spatial Reference Object for the input EPSG

        """
        if isinstance(epsg, (str, unicode)):
            try:
                epsg = int(epsg.replace('EPSG:', ''))
            except:
                raise TypeError('EPSG must be a number. Got - {}'.format(epsg))

        if not isinstance(epsg, (int, long)):
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

    def getEPSGFromSRS(self, osr_srs, bOnlineLookup=False, bUpdateToCorrectDefn=False):
        """Get the EPSG number for a Spatial Reference System.

           If the Spatial reference system does not contain an EPSG number it will attempt to find one.

           If bOnlineLookup is set to True it will use the online lookup service at  http://prj2epsg.org to attempt to
           identify the EPSG. This service will return a list of potential EPSG numbers. If a direct match is not found,
           then it will return None.

           If bOnlineLookup is set to False OR the online lookup fails, then it will attempt to match the input
           spatial reference system by looping through the list and attempting a proj4 string match if this still fails,
           then it will attempt a match to australian projected coordinate system

           adapted from : https://gis.stackexchange.com/a/8888

        Args:
            osr_srs (osgeo.osr.SpatialReference):
            bOnlineLookup (bool): Use the Open Geo Lookup Service
            bUpdateToCorrectDefn (bool): overwrite existing with the found definition
        Returns:
            int: EPSG Number for the matched Spatial Reference System.

        """
        warnings.warn('bOnlineLookup is deprecated', PendingDeprecationWarning)

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

            connection =  sqlite3.connect(crs_database)
            cursor = connection.cursor()
            cs_df = pd.read_sql_query('Select  * from  alias_name', connection)
            del cursor
            connection.close()

            cs_df = cs_df[cs_df['alt_name'].str.contains(srsName.replace('GDA',''), case=False)]

            if len(cs_df) == 1:
                epsg = int(cs_df['code'])

            if len(cs_df):
                LOGGER.debug('Matched to ESRI alias {}'.format(srsName))

        if bUpdateToCorrectDefn and epsg > 0:
            self.getFromEPSG(int(epsg))

        return int(epsg)


def getCRSfromRasterFile(raster_file):
    """Create a CRS object from a raster file.

    Args:
        raster_file (str): The path and filename to a raster file

    Returns:
        pyprecag.crs.crs: An object representing the raster file's coordinate system.

    """

    gdalRaster = gdal.Open(os.path.normpath(raster_file))
    rast_crs = crs()
    rast_crs.getFromWKT(gdalRaster.GetProjectionRef())
    del gdalRaster
    return rast_crs


def getCoordTransformation(inSR, outSR):
    """Get the coordinate transformation between two input spatial references
    Args:
        inSR (osgeo.osr.SpatialReference): Input Spatial Reference System
        outSR (osgeo.osr.SpatialReference): Output Spatial Reference System
    Returns:
        osgeo.osr.CoordinateTransformation:  The Coordinate Transformation object
    """
    return osr.CoordinateTransformation(inSR, outSR)


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

    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        if latitude < 0.0:
            return 0
        else:
            return 1

    utm_crs = osr.SpatialReference()
    utm_crs.SetWellKnownGeogCS("WGS84")  # Set geographic coordinate system to handle latitude/longitude
    zone = get_utm_zone(longitude)
    utm_crs.SetUTM(get_utm_zone(longitude), is_northern(latitude))
    wgs84_crs = utm_crs.CloneGeogCS()  # Clone ONLY the geographic coordinate system

    return zone, utm_crs, wgs84_crs


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

    # Based On :https://stackoverflow.com/a/10239676
    from pyproj import Proj, transform

    # Coordinates need to be in wgs84 so project them
    if xy_epsg != 4326:
        inProj = Proj(init='epsg:{}'.format(xy_epsg))
        outProj = Proj(init='epsg:4326')
        longitude, latitude = transform(inProj, outProj, x_coord, y_coord)
    else:
        longitude, latitude = x_coord, y_coord

    utm_zone = int(1 + (longitude + 180.0) / 6.0)

    # Determines if given latitude is a northern for UTM        1 is northern, 0 is southern
    is_northern = int(latitude > 0.0)

    utm_crs = crs()
    utm_crs.srs = osr.SpatialReference()

    # Use Australian or New Zealand local systems otherwise use global wgs84 UTM zones.
    if (108.0 <= longitude <= 155.0) and (-45.0 <= latitude <= -10.0):
        # set to Australian GDA94 MGA Zone XX
        utm_crs.srs.ImportFromEPSG(int('283{}'.format(utm_zone)))
    elif (166.33 <= longitude <= 178.6) and (-47.4 <= latitude <= -34.0):
        # set to NZGD2000 / New Zealand Transverse Mercator 2000
        utm_crs.srs.ImportFromEPSG(2193)
    else:
        # Set to Global WGS84 Utm Zone.
        utm_crs.srs.SetWellKnownGeogCS("WGS84")
        utm_crs.srs.SetUTM(utm_zone, is_northern)
        utm_crs.srs.AutoIdentifyEPSG()

    utm_crs.set_epsg(utm_crs.srs.GetAuthorityCode(None))

    utm_crs.epsg_predicted = False
    utm_crs.crs_wkt = utm_crs.srs.ExportToWkt()
    utm_crs.proj4 = utm_crs.srs.ExportToProj4()

    if utm_crs.srs.IsProjected():
        return utm_crs
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

    if not isinstance(distance_metres, (int, long, float)):
        raise TypeError('distance_metres must be a floating number.')

    # get the required Spatial reference systems
    utm_crs = getProjectedCRSForXY(longitude, latitude)
    wgs84SRS =  utm_crs.srs.CloneGeogCS()
    # zone, utmSRS, wgs84SRS = getUTMfromWGS84(longitude, latitude)

    # create transform component
    wgs842utm_transform = osr.CoordinateTransformation(wgs84SRS, utm_crs.srs)  # (<from>, <to>)
    # Project to UTM
    easting, northing, alt = wgs842utm_transform.TransformPoint(longitude, latitude, 0)

    # create transform component
    utm2wgs84_transform = osr.CoordinateTransformation(utm_crs.srs, wgs84SRS)  # (<from>, <to>)

    # Project back to WGS84 after adding distance to eastings. returns easting, northing, altitude
    newLong, newLat, alt = utm2wgs84_transform.TransformPoint(easting + distance_metres, northing, 0)

    # Create a point objects.
    origPoint = geometry.Point(longitude, latitude)
    adjPoint = geometry.Point(newLong, newLat)

    # Calculate Distance
    distance_dd = origPoint.distance(adjPoint)

    # if the input was negative make sure the output is too
    if distance_metres < 0:
        distance_dd = -distance_dd

    return distance_dd




