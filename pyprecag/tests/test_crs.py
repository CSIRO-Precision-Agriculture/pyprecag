import warnings
from unittest import TestCase
import os
import time
import shutil
import tempfile
import rasterio
import logging

from pyprecag.tests import make_dummy_data
from pyprecag.crs import crs, getProjectedCRSForXY, getCRSfromRasterFile, getUTMfromWGS84, \
    distance_metres_to_dd
from fiona.crs import from_epsg

this_dir = os.path.abspath(os.path.dirname(__file__))

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

EPSG_28354_WKT = ('PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia'
                  '_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],'
                  'TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,'
                  'AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG",'
                  '"9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],'
                  'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],'
                  'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
                  'PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
                  'AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

ESRI_54_WKT_1 = ('PROJCS["GDA_1994_MGA_Zone_54",GEOGCS["GCS_GDA_1994",'
                 'DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS_1980",6378137.0,'
                 '298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
                 'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
                 'PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",141.0],'
                 'PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],'
                 'UNIT["Meter",1.0]]')

ESRI_54_WKT_2 = ('PROJCS["GDA94_MGA_zone_54",GEOGCS["GCS_GDA_1994",DATUM["Geocentric_Datum_of'
                 '_Australia_1994",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich"' 
                 ',0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],'
                 'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],'
                 'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
                 'PARAMETER["false_northing",10000000],UNIT["Meter",1]]')

ESRI_54_SUTM_WKT = ('PROJCS["UTM Zone 54, Southern Hemisphere",GEOGCS["WGS 84",DATUM["WGS_1984",'
                    'SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],'
                    'AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
                    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
                    'AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],'
                    'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],'
                    'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
                    'PARAMETER["false_northing",10000000],UNIT["Meter",1]]')

NZ_WKT = ('PROJCS["NZGD2000 / New Zealand Transverse Mercator 2000",GEOGCS["NZGD2000",'
          'DATUM["New_Zealand_Geodetic_Datum_2000",SPHEROID["GRS 1980",6378137,298.257222101,'
          'AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6167"]],'
          'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,'
          'AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4167"]],PROJECTION["Transverse_Mercator"],'
          'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",173],'
          'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",1600000],'
          'PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
          'AUTHORITY["EPSG","2193"]]')


class TestCrsClass(TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestCrsClass, cls).setUpClass()

        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)
        cls.singletif, cls.multitif = make_dummy_data.make_dummy_tif_files(TmpDir)
        global testFailed
        testFailed = False

    @classmethod
    def tearDownClass(cls):
        if not testFailed:
            print ('Tests Passed .. Deleting {}'.format(TmpDir))
            shutil.rmtree(TmpDir)

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_getFromEPSG(self):
        test = crs()
        test.getFromEPSG(28354)
        self.assertEqual(test.epsg_number, 28354)
        self.assertEqual(test.epsg, from_epsg(test.epsg_number))
        self.assertEqual(EPSG_28354_WKT, test.crs_wkt)

        self.assertEqual('+proj=utm +zone=54 +south +ellps=GRS80 '
                         '+towgs84=0,0,0,0,0,0,0 +units=m +no_defs ',
                         test.proj4)

    def test_getFromWKT_GDA1(self):
        test = crs()
        test.getFromWKT(EPSG_28354_WKT)

        self.assertEqual(test.epsg_number, 28354)
        self.assertEqual(test.epsg, from_epsg(test.epsg_number))
        self.assertEqual(test.proj4, '+proj=utm +zone=54 +south +ellps=GRS80 '
                                     '+towgs84=0,0,0,0,0,0,0 +units=m +no_defs ')
        self.assertEqual(EPSG_28354_WKT, test.crs_wkt)

        test.getFromWKT(ESRI_54_WKT_2)
        self.assertEqual(28354, test.epsg_number)
        self.assertEqual(from_epsg(test.epsg_number), test.epsg)
        self.assertEqual('+proj=utm +zone=54 +south +ellps=GRS80 +units=m +no_defs ', test.proj4)
        self.assertEqual(ESRI_54_WKT_2, test.crs_wkt)

    def test_getUTMfromWGS84(self):
        result = getUTMfromWGS84(138.679870, -34.037740)
        self.assertEqual(54, result[0])
        self.assertEqual(ESRI_54_SUTM_WKT, result[1].ExportToWkt())

    def test_distance_metres_to_dd(self):
        dist = distance_metres_to_dd(138.822027994089, -34.4842175261199, 500)
        self.assertEqual(round(dist, 8), 0.00544233)

    def testGetFromWKTUsingLookup_GDA(self):

        with warnings.catch_warnings(record=True) as w:
            test = crs()
            test.getFromWKT(EPSG_28354_WKT)

            if len(w) == 1:
                print ('Timeout Occurred')
                self.assertIn('could not be reached. Timeout after 10 seconds', str(w[-1].message))
            else:
                self.assertEqual(EPSG_28354_WKT, test.crs_wkt)

                self.assertEqual(test.epsg, from_epsg(test.epsg_number))
                # self.assertEqual(test.proj4, '+proj=utm +zone=54 +south +ellps=GRS80'
                #                              ' +units=m +no_defs ')
                self.assertEqual(test.epsg_predicted, False)

            test = crs()
            test.getFromWKT(ESRI_54_WKT_1)

            if len(w) == 1:
                print ('Timeout Occurred')
                self.assertIn('could not be reached. Timeout after 10 seconds', str(w[-1].message))
            else:
                self.assertEqual(ESRI_54_WKT_1, test.crs_wkt)

                self.assertEqual(test.epsg, from_epsg(test.epsg_number))
                # self.assertEqual(test.proj4, '+proj=utm +zone=54 +south +ellps=GRS80'
                #                              ' +units=m +no_defs ')
                self.assertEqual(test.epsg_predicted, True)


            test = crs()
            test.getFromWKT(ESRI_54_WKT_2)

            if len(w) == 1:
                print ('Timeout Occurred')
                self.assertIn('could not be reached. Timeout after 10 seconds', str(w[-1].message))
            else:
                self.assertEqual(ESRI_54_WKT_2, test.crs_wkt)

                self.assertEqual(test.epsg, from_epsg(test.epsg_number))
                # self.assertEqual(test.proj4, '+proj=utm +zone=54 +south +ellps=GRS80'
                #                              ' +units=m +no_defs ')
                self.assertEqual(test.epsg_predicted, True)

    def test_getProjectedCRSForXY(self):
        # wgs84 edge of z54-55
        result = getProjectedCRSForXY(143.95231, -37.79412, 4326)
        self.assertEqual(EPSG_28354_WKT, result.srs.ExportToWkt())

        # Same point as above but in AGD66
        result = getProjectedCRSForXY(143.95099, -37.79561, 4202)
        self.assertEqual(result.epsg_number, 28354)

        # a New Zealand point
        result = getProjectedCRSForXY(169.796934, -44.380541)
        self.assertEqual(result.epsg_number, 2193)
        self.assertEqual(NZ_WKT, result.crs_wkt)

        # A northern hemisphere point
        result = getProjectedCRSForXY(143.95099, 37.79561, 4326)
        self.assertEqual(result.epsg_number, 32654)
        self.assertIn('PROJCS["UTM Zone 54, Northern Hemisphere",GEOGCS["WGS 84",'
                      'DATUM["WGS_1984",SPHEROID["WGS 84"',
                      result.crs_wkt)

        # a southern hemisphere point outside australia
        result = getProjectedCRSForXY(165.95099, -37.79561, 4326)
        self.assertEqual(result.epsg_number, 32758)

    def test_getRasterFileCrs(self):
        rast_crs = getCRSfromRasterFile(os.path.realpath(
            this_dir + r"/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif"))
        self.assertEqual(rast_crs.epsg_number, None)
        self.assertEqual(rast_crs.crs_wkt, None)

        rast_crs = getCRSfromRasterFile(os.path.normpath(self.singletif))
        self.assertEqual(rast_crs.epsg_number, 28354)
        self.assertEqual(rast_crs.crs_wkt, rasterio.crs.CRS.from_epsg(28354).wkt)
