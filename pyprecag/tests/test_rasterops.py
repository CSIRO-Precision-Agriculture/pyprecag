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
# TmpDir = r'C:\data\temp'
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

class test_rasterOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_rasterOps, cls).setUpClass()
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

    def test_RasterSnapExtent(self):
        xmin, ymin = 350127.023547, 6059756.84652457
        xmax, ymax = xmin + 3.142 * 200, ymin + 3.142 * 550

        # print(xmin, ymin, xmax, ymax)

        self.assertItemsEqual([350125.0, 6059755.0, 350760.0, 6061485.0], raster_ops.RasterSnapExtent(xmin, ymin, xmax, ymax, 5))
        self.assertItemsEqual([350127.0, 6059756.5, 350755.5, 6061485.0], raster_ops.RasterSnapExtent(xmin, ymin, xmax, ymax, 0.5))

    def test_rescale_singleband(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            rescaled = raster_ops.rescale(src, 0, 255)
            rescaled2 = raster_ops.rescale(src, 0, 5)
            outMeta = src.meta.copy()

        self.assertEquals(np.nanmin(rescaled), 0)
        self.assertEquals(np.nanmax(rescaled), 255)

        self.assertEquals(np.nanmin(rescaled2), 0)
        self.assertEquals(np.nanmax(rescaled2), 5)

        # with rasterio.open(os.path.normpath(os.path.join(out_dir, r'test_singleband_94mga54-rescale0-255.tif'), 'w', **outMeta) as out:
        #     out.write_band(1, rescaled)
        #
        # with rasterio.open(os.path.normpath(os.path.join(out_dir, r'test_singleband_94mga54-rescale0-5.tif'), 'w', **outMeta) as out:
        #     out.write_band(1, rescaled2)

    def test_rescale_multiband(self):
        with rasterio.open(os.path.normpath(self.multitif)) as src:
            rescaled = raster_ops.rescale(src, 0, 255, band_num=3)
            rescaled2 = raster_ops.rescale(src, 0, 5, band_num=3)

        self.assertEquals(int(np.nanmin(rescaled)), 0)
        self.assertEquals(int(np.nanmax(rescaled)), 255)

        self.assertEquals(int(np.nanmin(rescaled2)), 0)
        self.assertEquals(int(np.nanmax(rescaled2)), 5)

    def test_normalise(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            norm = raster_ops.normalise(src)
            outMeta = src.meta.copy()

        self.assertEquals(float(np.nanmax(norm)), 3.3085429668426514)
        self.assertEquals(float(np.nanmin(norm)), -2.758878707885742)
        # with rasterio.open(os.path.normpath(os.path.join(out_dir, r'test_singleband_94mga54-normalise.tif'), 'w', **outMeta) as out:
        #     out.write_band(1, norm)

    def test_SingleBand_focalstatistics(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanmean)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), 0.0363851934671402)
            self.assertEquals(float(np.nanmax(arr)), 6.991053581237793)
            self.assertEquals(col,  'mean5x5_test_singleband_94mga54')

            arr, col = raster_ops.focal_statistics(src, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), 2.632201585583971e-06)
            self.assertEquals(float(np.nanmax(arr)), 0.5194834470748901)
            self.assertEquals(col, 'cv3x3_test_singleband_94mga54')

            arr, col = raster_ops.focal_statistics(src, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), 5.0)
            self.assertEquals(float(np.nanmax(arr)), 81.0)
            self.assertEquals(col, 'pixelcount9x9_test_singleband_94mga54')

    def test_MultiBand_focalstatistics(self):
        with rasterio.open(os.path.normpath(self.multitif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanstd)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), 1.1673586413962767e-05)
            self.assertEquals(float(np.nanmax(arr)), 0.3216133415699005)
            self.assertEquals(col,  'std5x5bd1_test_3band_94mga54')

            arr, col = raster_ops.focal_statistics(src, band_num=3, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), -34.32146072387695)
            self.assertEquals(float(np.nanmax(arr)), 517.6610717773438)
            self.assertEquals(col, 'cv3x3bd3_test_3band_94mga54')

            arr, col = raster_ops.focal_statistics(src, band_num=2, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEquals(float(np.nanmin(arr)), 5.0)
            self.assertEquals(float(np.nanmax(arr)), 81.0)
            self.assertEquals(col, 'pixelcount9x9bd2_test_3band_94mga54')
