import shutil
import unittest
from geopandas import GeoDataFrame
from pyprecag.convert import convert_csv_to_points, numeric_pixelsize_to_string
import time
import os
import tempfile
import logging

pyFile = os.path.basename(__file__)

TmpDir = tempfile.gettempdir()
# TmpDir = r'C:\data\temp'
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

class test_convert(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_convert, cls).setUpClass()
        if not os.path.exists(TmpDir): os.mkdir(TmpDir)

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
        print("%s: %.3f" % (self.id(), t))

    def run(self, result=None):
        global testFailed
        unittest.TestCase.run(self, result) # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed=True

    def test_convert_csv_to_points_EastingNorthings(self):
        file = os.path.realpath(this_dir + "/data/area1_yield_file_ascii_wgs84.csv")
        epsg = 28354
        out_file = os.path.join(TmpDir, os.path.basename(file).replace('.csv', '_{}.shp'.format(epsg)))

        gdfData, gdfCRS = convert_csv_to_points(file, out_file, coord_columns=['Easting', 'Northing'],
                                                coord_columns_epsg=epsg)

        self.assertIsInstance(gdfData, GeoDataFrame)
        self.assertEqual(len(gdfData), 34309)
        self.assertEqual(gdfCRS.epsg_number, epsg)
        self.assertTrue(os.path.exists(out_file))

    def test_convert_csv_to_points_WGS84(self):
        file = os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv")
        out_epsg = 28354
        out_file = os.path.join(TmpDir, os.path.basename(file).replace('.csv', '_{}.shp'.format(out_epsg)))

        gdfData, gdfCRS = convert_csv_to_points(file, out_file, coord_columns_epsg=4326, out_epsg=out_epsg)

        self.assertTrue(os.path.exists(out_file))
        self.assertIsInstance(gdfData, GeoDataFrame)
        self.assertEqual(len(gdfData), 10000)
        self.assertEqual(list(set(gdfData.geom_type)), ['Point'])

        import numpy as np
        np.testing.assert_almost_equal(list(gdfData.total_bounds),
                               [598603.84418501, 6053251.95026914, 599822.96368475, 6054429.82552078], 4)

        self.assertEqual(gdfCRS.epsg_number, out_epsg)
        self.assertEqual(gdfCRS.crs_wkt,
                          'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

    def test_convert_csv_to_points_WGS84_GuessEPSG(self):
        file = os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv")

        out_file = os.path.join(TmpDir, os.path.basename(file).replace('.csv', '_guessepsg.shp'))
        gdfData, gdfCRS = convert_csv_to_points(file, out_file, coord_columns_epsg=4326, out_epsg=-1)
        self.assertTrue(os.path.exists(out_file))
        self.assertIsInstance(gdfData, GeoDataFrame)
        self.assertEqual(len(gdfData), 10000)
        self.assertEqual(list(set(gdfData.geom_type)), ['Point'])
        self.assertEqual(gdfCRS.epsg_number, 28354)
        self.assertEqual(gdfCRS.crs_wkt,
                          'PROJCS["GDA94 / MGA zone 54",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28354"]]')

    def test_numeric_pixelsize_to_string(self):
        self.assertEqual(numeric_pixelsize_to_string(0.42), '42cm')
        self.assertEqual(numeric_pixelsize_to_string(0.125), '125mm')
        self.assertEqual(numeric_pixelsize_to_string(2.0), '2m')
        self.assertEqual(numeric_pixelsize_to_string(6000),'6km')
        self.assertEqual(numeric_pixelsize_to_string(1500), '1500m')
