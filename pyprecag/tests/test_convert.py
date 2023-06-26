import shutil
import unittest

from pyprecag.tests import make_dummy_tif_files, setup_folder, KEEP_TEST_OUTPUTS
from pyprecag.tests.test_crs import EPSG_28354_WKT
from pyprecag.convert import *
from pyprecag import crs as pyprecag_crs

import time
import os
import tempfile
import logging
from shapely.geometry import LineString

import geopandas as gpd
import pandas as pd
import numpy.testing as npt

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
            for ea in ['pandas', 'numpy','gdal','geopandas','pyproj']:
                if 'gdal' == ea:
                    from osgeo import gdal

                else:
                    exec('import {}'.format( ea))

                if ea == 'pyproj':
                    exec('pyproj.show_versions()')
                else:
                    print('{:15}\t{}'.format(ea, eval('{}.__version__'.format(ea))))


            print('Tests Failed .. {}'.format(cls.TmpDir))

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = ''

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)


    def test_convert_csv_to_points_EastingNorthings(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area1_yield_ascii_wgs84.csv"))
        epsg = 28354
        self.test_outdir  = setup_folder(self.TmpDir, new_folder=self._testMethodName)
        out_file = os.path.join(self.test_outdir , os.path.basename(in_file).replace('.csv', '_{}.shp'.format(epsg)))

        gdf_data, gdf_crs = convert_csv_to_points(in_file, out_file,
                                                  coord_columns=['Easting', 'Northing'],
                                                  coord_columns_epsg=epsg)

        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(14756, len(gdf_data))
        self.assertEqual(epsg, gdf_crs.epsg_number)
        self.assertTrue(os.path.exists(out_file))

    def test_convert_csv_to_points_WGS84(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv"))
        out_epsg = 28354
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)
        out_file = os.path.join(self.test_outdir , os.path.basename(in_file).replace('.csv', '_{}.shp'.format(out_epsg)))

        gdf_data, gdf_crs = convert_csv_to_points(in_file, out_file, coord_columns_epsg=4326,
                                                  out_epsg=out_epsg)

        self.assertTrue(os.path.exists(out_file))
        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(1543, len(gdf_data))
        self.assertEqual(['Point'], list(set(gdf_data.geom_type)))

        import numpy as np
        np.testing.assert_almost_equal([598868.6709, 6054050.0529, 599242.9049, 6054415.3845],
                                       list(gdf_data.total_bounds), 4)

        self.assertEqual(out_epsg, gdf_crs.epsg_number)
        self.assertEqual(EPSG_28354_WKT[:154], gdf_crs.crs_wkt[:154])

    def test_convert_csv_to_points_WGS84_GuessEPSG(self):
        in_file = os.path.realpath(os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv"))

        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)
        self.test_outdir = os.path.join(out_fold, os.path.basename(in_file).replace('.csv', '_guessepsg.shp'))

        gdf_data, gdf_crs = convert_csv_to_points(in_file, self.test_outdir , coord_columns_epsg=4326, out_epsg=-1)

        self.assertTrue(os.path.exists(self.test_outdir))
        self.assertIsInstance(gdf_data, gpd.GeoDataFrame)
        self.assertEqual(1543, len(gdf_data))
        self.assertEqual(['Point'], list(set(gdf_data.geom_type)) )
        self.assertEqual(28354, gdf_crs.epsg_number)
        self.assertEqual(EPSG_28354_WKT[:154], gdf_crs.crs_wkt[:154])

    def test_numeric_pixelsize_to_string(self):
        self.assertEqual('42cm', numeric_pixelsize_to_string(0.42))
        self.assertEqual('125mm', numeric_pixelsize_to_string(0.125))
        self.assertEqual('2m', numeric_pixelsize_to_string(2.0))
        self.assertEqual('6km', numeric_pixelsize_to_string(6000))
        self.assertEqual('1500m', numeric_pixelsize_to_string(1500))
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

        c_line_gdf['geom_noZ'] = c_line_gdf['geometry'].apply(lambda x: drop_z(x))

        self.assertEqual(c_line_gdf['geometry'].iloc[1].wkt, c_line_gdf['geom_noZ'].iloc[1].wkt)
        self.assertNotEqual(c_line_gdf['geometry'].iloc[0].wkt, c_line_gdf['geom_noZ'].iloc[0].wkt)

        self.assertFalse(c_line_gdf['geom_noZ'].iloc[0].has_z)
        self.assertFalse(c_line_gdf['geom_noZ'].iloc[-1].has_z)
