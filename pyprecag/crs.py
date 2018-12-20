import json
import logging
import math
import os
import warnings
from urllib import urlencode, urlopen

from fiona.crs import from_string, from_epsg
from osgeo import osr, gdal
from shapely import geometry

from . import config
from .errors import SpatialReferenceError

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
# LOGGER.setLevel("DEBUG")
DEBUG = config.get_config_key('debug_mode')  # LOGGER.isEnabledFor(logging.DEBUG)


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
            SRS attributes to match the official EPGS database attributes
        Args:
            crs_wkt (unicode):    A string representing the coordinate system well known text
            bUpdateByEPSG (bool): If True and the epsg_number is predicted, overwrite crs attributes with those from
                                  official epsg_number database.
        Examples:
            >>> coordSys = crs()
            >>> coordSys.getFromWKT('PROJCS["GDA_1994_MGA_Zone_54",GEOGCS["GCS_GDA_1994",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",141.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
            >>> coordSys.crs_wkt
            'PROJCS["GDA_1994_MGA_Zone_54",GEOGCS["GCS_GDA_1994",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",141.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
            >>> coordSys.epsg_number
            28354
            >>> coordSys.epsg_predicted
            False
            >>> coordSys.proj4
            '+proj=utm +zone=54 +south +ellps=GRS80 +units=m +no_defs '

        """
        if crs_wkt == '':
            warnings.warn('crs_wkt is blank. No coordinate information set')
            # LOGGER.warning('WARNING: crs_wkt is blank. No Coordinate Information Set')
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
                code = self.getEPSGFromSRS(source, bOnlineLookup=True, bUpdateToCorrectDefn=bUpdateByEPSG)
            else:
                code = int(source.GetAuthorityCode(None))

            if code is not None:
                self.set_epsg(code)


    def getFromEPSG(self, epsg):
        # noinspection PyPep8
        """Create OGR Spatial Reference Object for an epsg_number number
                        Args:
                            epsg (int):    A Valid EPSG Number

                        Returns:
                            osgeo.osr.SpatialReference: The Spatial Reference Object for the input EPSG

                        Examples:
                            >>> coordSys = crs()
                            >>> coordSys.getFromEPSG(28354)
                            >>> coordSys.crs_wkt
                            'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]'
                            >>> coordSys.proj4
                            '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '

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

    # noinspection SpellCheckingInspection
    def getEPSGFromSRS(self, osr_srs, bOnlineLookup=False, bUpdateToCorrectDefn=False):
        """Get the EPSG number for a Spatial Reference System.

           If the Spatial reference system does not contain an EPSG number it will attempt to find one.

           If bOnlineLookup is set to True it will use the online lookup service at  http://prj2epsg.org to attempt to
           identify the EPSG. This service will return a list of potential EPSG numbers. If a direct match is not found,
           then it will return None.

           If bOnlineLookup is set to False OR the online lookup fails, then it will attempt to match the input
           spatial reference system by looping through the list and attempting a proj4 string match if this still fails,
           then it will attempt a match to australian projected coordinate system

            # adapated from : https://gis.stackexchange.com/a/8888

        Args:
            osr_srs (osgeo.osr.SpatialReference):
            bOnlineLookup (bool): Use the Open Geo Lookup Service
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

        # if osr_srs.IsProjected():
        #     epsg_number = osr_srs.GetAttrValue("PROJCS|AUTHORITY", 1)
        # else:
        #     epsg_number = osr_srs.GetAttrValue("PRIMEM|AUTHORITY", 1)

        # Then through a lookup
        if epsg is None and bOnlineLookup:
            self.epsg_predicted = True
            query = urlencode({
                'exact': True,
                'error': True,
                'mode': 'wkt',
                'terms': osr_srs.ExportToWkt()})
            try:
                webres = None
                crsURL = config.get_config_key('crsLookupURL')
                LOGGER.debug('Checking against OpenGeo service ({})'.format(crsURL))
                webres = urlopen(crsURL, query)
            except:
                LOGGER.warning('WARNING: OpenGeo service ({}) could not be reached. '.format(crsURL))

            if webres is not None:
                jres = json.loads(webres.read())
                if len(jres['codes']) == 1:
                    epsg = jres['codes'][0]['code']
                    LOGGER.debug('Matched to EPSG {} {}'.format(jres['codes'][0]['code'], jres['codes'][0]['name']))
                    self.epsg_predicted = False
                elif len(jres['codes']) > 1:
                    LOGGER.debug(
                        '\n\nEPSG lookup found {} matches. Attempting to refine by comparing proj4 strings'.format(
                            len(jres['codes'])))
                    for i in reversed(range(len(jres['codes']))):
                        epsg = None
                        tmpSrs = osr.SpatialReference()
                        res = tmpSrs.ImportFromEPSG(int(jres['codes'][i]['code']))

                        if res != 0:
                            raise RuntimeError(repr(res) + ': could not import from EPSG')

                        # create a dictionary mapping using fiona.crs.from_string to ensure elements are in
                        # the same order.
                        tmpProj4Dict = from_string(tmpSrs.ExportToProj4())

                        if from_string(osr_srs.ExportToProj4()) == tmpProj4Dict:
                            epsg = jres['codes'][i]['code']
                        else:
                            # remove towgs84 value if all 0's as it is not always implemented yet for gda2020
                            if 'towgs84' in tmpProj4Dict:
                                if tmpProj4Dict['towgs84'] == '0,0,0,0,0,0,0':
                                    del tmpProj4Dict['towgs84']

                            if from_string(osr_srs.ExportToProj4()) == tmpProj4Dict:
                                epsg = jres['codes'][i]['code']

                        if epsg is None:
                            del jres['codes'][i]

                    if len(jres['codes']) == 1:
                        epsg = jres['codes'][0]['code']
                        LOGGER.debug('Refined match returns EPSG {} {}'.format(jres['codes'][0]['code'],
                                                                               jres['codes'][0]['name']))
                    else:
                        mess = 'ERROR:-EPSG lookup found {} matches. Please properly define the projection ' \
                               'and try again.\nThe matches were:'.format(len(jres['codes']))
                        for i in range(len(jres['codes'])):
                            mess = mess + '\t{0:>7}    {1}'.format(jres['codes'][i]['code'], jres['codes'][i]['name'])
                        raise SpatialReferenceError(mess)

        # if still none, then attempt to map for common aussie prj's only
        if epsg is None or epsg == 0:
            self.epsg_predicted = True
            if osr_srs.IsProjected():
                srsName = osr_srs.GetAttrValue("PROJCS", 0)
            else:
                srsName = osr_srs.GetAttrValue("GEOGCS", 0)

            LOGGER.debug('Attempting to map {} to Australian GDA projections and WGS84.'.format(srsName), )
            # TODO: Implement GDA2020
            if osr_srs.IsProjected() and 'GDA' in srsName.upper() and 'MGA' in srsName.upper():
                epsg = 28300 + int(srsName[-2:])
            elif srsName == 'GCS_WGS_1984':
                epsg = 4326
            elif srsName == 'GCS_GDA_1994':
                epsg = 4283
            if epsg is None:
                epsg = 0
                LOGGER.warning('WARNING: No EPSG match found')
            else:
                LOGGER.debug('Aussie match found: EPSG:{}'.format(epsg))

        if bUpdateToCorrectDefn and epsg > 0:
            self.getFromEPSG(int(epsg))

        return int(epsg)


def getCRSfromRasterFile(raster_file):
    """Create a CRS object from a raster file.

    Args:
        raster_file (str): The path and filename to a raster file

    Returns:
        pyprecag.crs.crs: An object representing the raster file's coordinate system.

    Example:
        >>> rast_crs = getCRSfromRasterFile(r'../test/data/test_singleband_94mga54.tif')
        >>> rast_crs.epsg_number
        28354
        >>> rast_crs.epsg_predicted
        False
        >>> rast_crs.crs_wkt
        'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]'

    """

    gdalRaster = gdal.Open(os.path.normpath(raster_file))
    rast_crs = crs()
    rast_crs.getFromWKT(gdalRaster.GetProjectionRef())
    del gdalRaster
    return rast_crs


# def getEPSGForLongitude(longitude):
#     """Given a longitude, return the EPSG code of the MGA/GDA94 Zone
#
#     Args:
#         longitude (float):  A Longitude value
#
#     Returns:
#         zone: the MGA Zone
#         int: EPSG Number
#     Example:
#         >>> getEPSGForLongitude(138.679870)
#         (54, 28354)
#     """
#
#     Zone = 50 + math.trunc((longitude - 114) / 6)
#     # if debugPrint: print('getLongGDASpatRef: Derived UTM zone: {}'.format(str(Zone))
#
#     # EPSG Codes for GDA UTM zones
#     # 28300 + Zone number. i.e. 28350 = GDA94 / MGA zone 50
#     EPSGCode = 28300 + Zone
#     # if debugPrint: print('GDA EPSG: {}'.format(EPSGCode))
#     return Zone, EPSGCode


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

    Examples:
        >>> result = getUTMfromWGS84(138.679870, -34.037740)
        >>> result[0]
        54
        >>> result[1].ExportToWkt()
        'PROJCS["UTM Zone 54, Southern Hemisphere",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["Meter",1]]'
    """

    # Add Function + arguments to Log
    # [logger.info(ea) for ea in general.print_functions_string(inspect.currentframe(), self.__class__.__name__)]
    # print('\n'.join(general.print_functions_string(inspect.currentframe())))

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

        Examples:
            >>> result = getProjectedCRSForXY(143.95231,-37.79412,4326)    # wgs84 edge of z54-55
            >>> result.srs.ExportToWkt()
            'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]'
            >>> result = getProjectedCRSForXY(143.95099,-37.79561,4202)    # Same point as above but in AGD66
            >>> result.crs_wkt
            'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]'

        """

    # Add Function + arguments to Log
    # [logger.info(ea) for ea in general.print_functions_string(inspect.currentframe(), self.__class__.__name__)]
    # print('\n'.join(general.print_functions_string(inspect.currentframe())))

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

    # if in australia use GDA94 MGA zones otherwise use UTM system
    if (111.927 < longitude < 154.047) and (-44.5362 < latitude < -8.69876):
        utm_crs.srs.ImportFromEPSG(int('283{}'.format(utm_zone)))
    else:
        utm_crs.srs.SetWellKnownGeogCS("WGS84")  # Set geographic coordinate system to handle latitude/longitude
        utm_crs.srs.SetUTM(utm_zone, is_northern)
        utm_crs.srs.AutoIdentifyEPSG()

    utm_crs.set_epsg(utm_crs.srs.GetAuthorityCode(None))

    # wgs84_crs = utm_crs.CloneGeogCS()  # Clone ONLY the geographic coordinate system

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

    Examples:
        >>> dist = distance_metres_to_dd (138.822027994089, -34.4842175261199, 500)
        >>> round(dist,8)
        0.00544233
    """

    # Add Function + arguments to Log
    # [logger.info(ea) for ea in general.print_functions_string(inspect.currentframe(), self.__class__.__name__)]
    # print('\n'.join(general.print_functions_string(inspect.currentframe())))

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

    # if the input was negative then keep output aggregateDist_m
    if distance_metres < 0:
        distance_dd = -distance_dd
    # print('Decimal Degrees Distance for {} m is {:f}'.format(distance_metres,distance_dd))

    return distance_dd

#
# def distanceBetween_wgs84_Points_in_metres(point1, point2):
#     """Calculate the distance in metres between two WGS84 points
#     It assumes the input points are in WGS84
#
#     Args:
#         point1 (shapely.geometry.Point):
#         point2 (shapely.geometry.Point):
#
#     Returns:
#         float: Distance in metres
#     """
#
#     # Add Function + arguments to Log
#     # [logger.info(ea) for ea in general.print_functions_string(inspect.currentframe(), self.__class__.__name__)]
#     # print('\n'.join(general.print_functions_string(inspect.currentframe())))
#
#     zone, utmSRS, wgs84SRS = getUTMfromWGS84(*point1['geometry']['coordinates'])
#
#     # create transform component
#     easting, northing, alt = osr.CoordinateTransformation(wgs84SRS, utmSRS).TransformPoint(
#         *point1['geometry']['coordinates'])
#     easting2, northing2, alt2 = osr.CoordinateTransformation(wgs84SRS, utmSRS).TransformPoint(
#         *point2['geometry']['coordinates'])
#
#     origPoint = geometry.Point(easting, northing)
#     adjPoint = geometry.Point(easting2, northing2)
#
#     distance_dd = origPoint.distance(adjPoint)
#     return distance_dd
