import os
import shutil
import tempfile
import time
import unittest

import pandas as pd

from pyprecag.convert import convert_csv_to_points
from pyprecag.tests import setup_folder, KEEP_TEST_OUTPUTS
from pyprecag.vector_ops import thin_point_by_distance

PY_FILE = os.path.basename(__file__)
TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data')


class test_vectorOps(unittest.TestCase):
    failedTests = []
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_vectorOps, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Folder Exists.. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

    def run(self, result=None):
        
        unittest.TestCase.run(self, result) # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            self.failedTests=True

    def test_thinPointByDistance_wgs84(self):
        file = os.path.realpath(os.path.join(THIS_DIR, "area1_yield_ascii_wgs84.csv"))
        out_epsg = 28354
        out_file = os.path.join(self.TmpDir, os.path.basename(file).replace('.csv', '_{}.shp'.format(out_epsg)))
        gdf, gdfCRS = convert_csv_to_points(file, out_file, coord_columns_epsg=4326, out_epsg=out_epsg)

        result = thin_point_by_distance(gdf, gdfCRS, 2.5)
        stats = result.groupby(by='filter', dropna=False).agg(count=pd.NamedAgg(column='filter', aggfunc='count'))

        pd.testing.assert_frame_equal(pd.DataFrame.from_records([{'filter': 'pointXY (2.5m)', 'count': 1338},
                                                                 {'filter': pd.NA, 'count': 0}], index='filter'), stats)

        self.assertEqual(len(result), len(gdf))
        self.assertEqual(result.crs, gdf.crs)


    def test_thinPointByDistance_mga54(self):

        file = os.path.realpath(os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv"))
        out_epsg = 28354
        out_file = os.path.join(self.TmpDir, os.path.basename(file).replace('.csv', '_{}.shp'.format(out_epsg)))
        gdf, gdfCRS = convert_csv_to_points(file, out_file, coord_columns_epsg=4326, out_epsg=out_epsg)

        result = thin_point_by_distance(gdf, gdfCRS, 2.5)
        stats = result.groupby(by='filter', dropna=False).agg(count=pd.NamedAgg(column='filter', aggfunc='count'))
        pd.testing.assert_frame_equal(pd.DataFrame.from_records([{'filter': 'pointXY (2.5m)', 'count': 109},
                                                                 {'filter': pd.NA, 'count': 0}], index='filter'), stats)

        self.assertEqual(len(result), len(gdf))
        self.assertEqual(result.crs, gdf.crs)
