import os
import shutil
import tempfile
import time
import unittest

from pyprecag.convert import convertCsvToPoints
from pyprecag.vector_ops import thin_point_by_distance

pyFile = os.path.basename(__file__)

TmpDir = tempfile.gettempdir()
# TmpDir = r'C:\data\temp'
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

class test_vectorOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_vectorOps, cls).setUpClass()
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

    def test_thin_point_by_distance_mga54(self):

        file = os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv")
        out_epsg = 28354
        out_file = os.path.join(TmpDir, os.path.basename(file).replace('.csv', '_{}.shp'.format(out_epsg)))
        gdf, gdfCRS = convertCsvToPoints(file, out_file, coord_columns_EPSG=4326, out_EPSG=out_epsg)

        result = thin_point_by_distance(gdf, gdfCRS, 2.5)
        self.assertEquals(len(result), len(gdf))
        self.assertEquals(len(result[result['filter'].isnull()]), 8406)
        self.assertEquals(result.crs, gdf.crs)
