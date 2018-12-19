import os
import shutil
import tempfile
import time
import unittest
import rasterio

from pyprecag.tests import make_dummy_data
from pyprecag import raster_ops
import numpy as np

pyFile = os.path.basename(__file__)

TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])


class test_rasterOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_rasterOps, cls).setUpClass()
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

    def test_rasterSnapExtent(self):
        xmin, ymin = 350127.023547, 6059756.84652457
        xmax, ymax = xmin + 3.142 * 200, ymin + 3.142 * 550

        self.assertItemsEqual([350125.0, 6059755.0, 350760.0, 6061485.0], raster_ops.raster_snap_extent(xmin, ymin, xmax, ymax, 5))
        self.assertItemsEqual([350127.0, 6059756.5, 350755.5, 6061485.0], raster_ops.raster_snap_extent(xmin, ymin, xmax, ymax, 0.5))

    def test_rescaleSingleBand(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            rescaled = raster_ops.rescale(src, 0, 255)
            rescaled2 = raster_ops.rescale(src, 0, 5)
            outMeta = src.meta.copy()

        self.assertEqual(np.nanmin(rescaled), 0)
        self.assertEqual(np.nanmax(rescaled), 255)

        self.assertEqual(np.nanmin(rescaled2), 0)
        self.assertEqual(np.nanmax(rescaled2), 5)

    def test_rescaleMultiBand(self):
        with rasterio.open(os.path.normpath(self.multitif)) as src:
            rescaled = raster_ops.rescale(src, 0, 255, band_num=3)
            rescaled2 = raster_ops.rescale(src, 0, 5, band_num=3)

        self.assertEqual(int(np.nanmin(rescaled)), 0)
        self.assertEqual(int(np.nanmax(rescaled)), 255)

        self.assertEqual(int(np.nanmin(rescaled2)), 0)
        self.assertEqual(int(np.nanmax(rescaled2)), 5)

    def test_normalise(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            norm = raster_ops.normalise(src)

        self.assertEqual(float(np.nanmax(norm)), 3.3085429668426514)
        self.assertEqual(float(np.nanmin(norm)), -2.758878707885742)

    def test_focalStatisticsSingleBand(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanmean)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), 0.0363851934671402)
            self.assertEqual(float(np.nanmax(arr)), 6.991053581237793)
            self.assertEqual(col,  'mean5x5_test_singleband_94mga54')

            arr, col = raster_ops.focal_statistics(src, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), 2.632201585583971e-06)
            self.assertEqual(float(np.nanmax(arr)), 0.5194834470748901)
            self.assertEqual(col, 'cv3x3_test_singleband_94mga54')

            arr, col = raster_ops.focal_statistics(src, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), 5.0)
            self.assertEqual(float(np.nanmax(arr)), 81.0)
            self.assertEqual(col, 'pixelcount9x9_test_singleband_94mga54')

    def test_focalStatisticsMultiBand(self):
        with rasterio.open(os.path.normpath(self.multitif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanstd)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), 1.1673586413962767e-05)
            self.assertEqual(float(np.nanmax(arr)), 0.3216133415699005)
            self.assertEqual(col,  'std5x5bd1_test_3band_94mga54')

            arr, col = raster_ops.focal_statistics(src, band_num=3, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), -34.32146072387695)
            self.assertEqual(float(np.nanmax(arr)), 517.6610717773438)
            self.assertEqual(col, 'cv3x3bd3_test_3band_94mga54')

            arr, col = raster_ops.focal_statistics(src, band_num=2, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEqual(float(np.nanmin(arr)), 5.0)
            self.assertEqual(float(np.nanmax(arr)), 81.0)
            self.assertEqual(col, 'pixelcount9x9bd2_test_3band_94mga54')
