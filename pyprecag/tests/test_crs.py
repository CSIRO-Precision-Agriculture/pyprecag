from unittest import TestCase
import os
import time
import shutil
import tempfile
import rasterio
from pyprecag.tests import make_dummy_data
from pyprecag.crs import crs, getProjectedCRSForXY, getCRSfromRasterFile
from fiona.crs import from_epsg


this_dir = os.path.abspath(os.path.dirname(__file__))

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

class TestCrsClass(TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestCrsClass, cls).setUpClass()
        if not os.path.exists(TmpDir): os.mkdir(TmpDir)
        cls.singletif, cls.multitif = make_dummy_data.make_dummy_tif_files(TmpDir)
        global testFailed
        testFailed = False

    @classmethod
    def tearDownClass(cls):
        if not testFailed:
            print 'Deleting folder {}'.format(TmpDir)
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
        self.assertEqual(test.crs_wkt,
                          'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')
        self.assertEqual(test.proj4,
                          '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ')

    def test_getFromWKT_GDA1(self):
        test = crs()
        test.getFromWKT(
            'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

        self.assertEqual(test.epsg_number, 28354)
        self.assertEqual(test.epsg, from_epsg(test.epsg_number))
        self.assertEqual(test.proj4,
                          '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ')
        self.assertEqual(test.crs_wkt,
                          'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

    def test_getFromWKT_GDA2(self):
        test = crs()
        test.getFromWKT(
            'PROJCS["GDA_1994_MGA_Zone_54",GEOGCS["GCS_GDA_1994",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",141.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
        self.assertEqual(test.crs_wkt,
                          'PROJCS["GDA_1994_MGA_Zone_54",GEOGCS["GCS_GDA_1994",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",141.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
        self.assertEqual(test.epsg, from_epsg(test.epsg_number))
        self.assertEqual(test.proj4, '+proj=utm +zone=54 +south +ellps=GRS80 +units=m +no_defs ')

    def test_getProjectedCRSForXY(self):
        result = getProjectedCRSForXY(143.95231, -37.79412, 4326)  # wgs84 edge of z54-55
        self.assertEqual(result.srs.ExportToWkt(),
                          'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

        result = getProjectedCRSForXY(143.95099, -37.79561, 4202)  # Same point as above but in AGD66
        self.assertEqual(result.epsg_number, 28354)

        # A northern hemisphere point
        result = getProjectedCRSForXY(143.95099, 37.79561, 4202)
        self.assertEqual(result.epsg_number, 32654)
        self.assertEqual(result.crs_wkt, 'PROJCS["UTM Zone 54, Northern Hemisphere",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1],AUTHORITY["EPSG","32654"]]')

        # a southern hemisphere point outside australia
        result = getProjectedCRSForXY(165.95099, -37.79561, 4202)
        self.assertEqual(result.epsg_number, 32758)


    def test_getRasterFileCrs(self):

        rast_crs = getCRSfromRasterFile(os.path.realpath(this_dir + r"/data/area1_rgbi_jan_50cm_84sutm54.tif"))
        self.assertEqual(rast_crs.epsg_number,None)
        self.assertEqual(rast_crs.crs_wkt, None)
        #'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

        rast_crs = getCRSfromRasterFile(os.path.normpath(self.singletif))
        self.assertEqual(rast_crs.epsg_number,28354)
        self.assertEqual(rast_crs.crs_wkt, rasterio.crs.CRS.from_epsg(28354).wkt)
