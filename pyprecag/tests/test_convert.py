import shutil
import unittest

from pyprecag.tests import make_dummy_tif_files, setup_folder, KEEP_TEST_OUTPUTS

from pyprecag.convert import *
from pyprecag import crs as pyprecag_crs
import time
import os
import tempfile
import logging
from shapely.geometry import LineString
import geopandas as gpd
from pyprecag.tests.test_crs import EPSG_28354_WKT

PY_FILE = os.path.basename(__file__)
TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")


class TestConvert(unittest.TestCase):
    failedTests = []
    
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestConvert, cls).setUpClass()
        
        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)
        else:
            print('Tests Failed .. {}'.format(cls.TmpDir))

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

    def run(self, result=None):

        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)


    def test_convert_csv_to_points_EastingNorthings(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area1_yield_ascii_wgs84.csv"))
        epsg = 28354
        out_file = os.path.join(self.TmpDir, os.path.basename(in_file).replace('.csv', '_{}.shp'.format(epsg)))

        gdf_data, gdf_crs = convert_csv_to_points(in_file, out_file,
                                                  coord_columns=['Easting', 'Northing'],
                                                  coord_columns_epsg=epsg)

        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(len(gdf_data), 13756)
        self.assertEqual(gdf_crs.epsg_number, epsg)
        self.assertTrue(os.path.exists(out_file))

    def test_convert_csv_to_points_WGS84(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv"))
        out_epsg = 28354
        out_file = os.path.join(self.TmpDir, os.path.basename(in_file).replace(
            '.csv', '_{}.shp'.format(out_epsg)))

        gdf_data, gdf_crs = convert_csv_to_points(in_file, out_file, coord_columns_epsg=4326,
                                                  out_epsg=out_epsg)

        self.assertTrue(os.path.exists(out_file))
        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(len(gdf_data), 1543)
        self.assertEqual(list(set(gdf_data.geom_type)), ['Point'])

        import numpy as np
        np.testing.assert_almost_equal(list(gdf_data.total_bounds),
                                       [598868.6709, 6054050.0529, 599242.9049, 6054415.3845], 4)

        self.assertEqual(gdf_crs.epsg_number, out_epsg)
        self.assertEqual(gdf_crs.crs_wkt[:154], EPSG_28354_WKT[:154])

    def test_convert_csv_to_points_WGS84_GuessEPSG(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv"))

        out_file = os.path.join(self.TmpDir, os.path.basename(in_file).replace('.csv', '_guessepsg.shp'))
        gdf_data, gdf_crs = convert_csv_to_points(in_file, out_file, coord_columns_epsg=4326,
                                                  out_epsg=-1)
        self.assertTrue(os.path.exists(out_file))
        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(len(gdf_data), 1543)
        self.assertEqual(list(set(gdf_data.geom_type)), ['Point'])
        self.assertEqual(gdf_crs.epsg_number, 28354)
        self.assertEqual(gdf_crs.crs_wkt[:154], EPSG_28354_WKT[:154])

    def test_numeric_pixelsize_to_string(self):
        self.assertEqual(numeric_pixelsize_to_string(0.42), '42cm')
        self.assertEqual(numeric_pixelsize_to_string(0.125), '125mm')
        self.assertEqual(numeric_pixelsize_to_string(2.0), '2m')
        self.assertEqual(numeric_pixelsize_to_string(6000), '6km')
        self.assertEqual(numeric_pixelsize_to_string(1500), '1500m')

    def test_cardinal_direction(self):
        self.assertEqual('E', deg_to_8_compass_pts(89.9999))
        self.assertEqual('N', deg_to_8_compass_pts(-5))
        self.assertEqual('SW', deg_to_8_compass_pts(210))
        self.assertEqual('N', deg_to_8_compass_pts(0.05))

        self.assertEqual('E', deg_to_16_compass_pts(89.9999))
        self.assertEqual('N', deg_to_16_compass_pts(-5))
        self.assertEqual('SSW', deg_to_16_compass_pts(210))
        self.assertEqual('N', deg_to_16_compass_pts(0.05))
        self.assertEqual('NE', deg_to_16_compass_pts(44))
        self.assertEqual('NW', deg_to_16_compass_pts(-44))

    def test_pts_linebearing(self):
        # line_bearing uses point_to_point_bearing so this doubles as it's tests as well
        c_line_gdf = gpd.GeoDataFrame({'geometry': [LineString([(740873, 6169764), (741269, 6169764)]),
                                                    LineString([(741000, 6169800), (741003, 6170012)]),
                                                    LineString([(741401, 6169800), (741057, 6170012)]),
                                                    LineString([(740900, 6169912), (740979, 6170094)])],
                                       'LineID'  : [1, 2, 3, 4]}, crs=pyprecag_crs.from_epsg(28354))

        c_line_gdf = GeoDataFrame(c_line_gdf, crs=pyprecag_crs.from_epsg(28354))

        self.assertEqual(90.0, line_bearing(c_line_gdf.iloc[0]))
        self.assertAlmostEqual(0.8107352192623694, line_bearing(c_line_gdf.iloc[1]), 4)
        self.assertAlmostEqual(301.64465903709254, line_bearing(c_line_gdf.iloc[2]), 4)
        self.assertAlmostEqual(23.464022404707464, line_bearing(c_line_gdf.iloc[3]), 4)

    def test_drop_z(self):
        c_line_gdf = gpd.GeoDataFrame({'geometry': [LineString([(740873, 6169764, 5), (741269, 6169764, 5)]),
                                                    LineString([(741000, 6169800), (741003, 6170012)]),
                                                    LineString([(741401.415, 6169800, 4), (741057, 6170012, 3)]),
                                                    LineString([(740900.861, 6169912, 2), (740979, 6170094, 5)])],
                                       'LineID'  : [1, 2, 3, 4]}, crs=pyprecag_crs.from_epsg(28354))

        c_line_gdf['geometry2'] = c_line_gdf['geometry'].apply(lambda x: drop_z(x))

        self.assertEqual(c_line_gdf['geometry'].iloc[1], c_line_gdf['geometry2'].iloc[1])
        self.assertNotEqual(c_line_gdf['geometry'].iloc[0], c_line_gdf['geometry2'].iloc[0])

        self.assertFalse(c_line_gdf['geometry2'].iloc[0].has_z)
        self.assertFalse(c_line_gdf['geometry2'].iloc[-1].has_z)