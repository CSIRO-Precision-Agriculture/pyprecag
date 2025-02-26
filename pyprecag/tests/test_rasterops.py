import os
import shutil
import tempfile
import time
import unittest
import rasterio
import six

from pyprecag.tests import make_dummy_tif_files, setup_folder, KEEP_TEST_OUTPUTS
from pyprecag import raster_ops
import numpy as np

PY_FILE = os.path.basename(__file__)

TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])


class test_rasterOps(unittest.TestCase):
    failedTests = []
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_rasterOps, cls).setUpClass()
        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print ('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = ''

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failures or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_rasterSnapExtent(self):
        xmin, ymin = 350127.023547, 6059756.84652457
        xmax, ymax = xmin + 3.142 * 200, ymin + 3.142 * 550

        six.assertCountEqual(
            self, [350125.0, 6059755.0, 350760.0, 6061485.0],
            raster_ops.raster_snap_extent(xmin, ymin, xmax, ymax, 5)
        )
        six.assertCountEqual(
            self, [350127.0, 6059756.5, 350755.5, 6061485.0],
            raster_ops.raster_snap_extent(xmin, ymin, xmax, ymax, 0.5)
        )

    def test_rescaleSingleBand(self):
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        with rasterio.open(os.path.normpath(self.singletif)) as src:
            rescaled0_255 = raster_ops.rescale(src, 0, 255)
            rescaled0_5 = raster_ops.rescale(src, 0, 5)
            out_meta = src.meta.copy()

        if KEEP_TEST_OUTPUTS:
            out_rescale = os.path.join(self.test_outdir,  'output_singleband_rescale0-255.tif')
            with rasterio.open(os.path.normpath(out_rescale), 'w', **out_meta) as out:
                out.write_band(1, rescaled0_255)

            out_rescale = os.path.join(self.test_outdir, 'output_singleband_rescale0-5.tif')
            with rasterio.open(os.path.normpath(out_rescale), 'w', **out_meta) as out:
                out.write_band(1, rescaled0_5)

        self.assertEqual(0, np.nanmin(rescaled0_255), )
        self.assertEqual(255, np.nanmax(rescaled0_255),)

        self.assertEqual(0, np.nanmin(rescaled0_5), )
        self.assertEqual(5, np.nanmax(rescaled0_5), )

    def test_rescaleMultiBand(self):
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        with rasterio.open(os.path.normpath(self.multitif)) as src:
            rescaled0_255 = raster_ops.rescale(src, 0, 255, band_num=3)
            rescaled0_5 = raster_ops.rescale(src, 0, 5, band_num=3)
            out_meta = src.meta.copy()

        out_meta.update({'count': 1})
        if KEEP_TEST_OUTPUTS:
            out_rescale = os.path.join(self.test_outdir,  'output_multiband_rescale0-255.tif')
            with rasterio.open(os.path.normpath(out_rescale), 'w', **out_meta) as out:
                out.write_band(1, rescaled0_255)

            out_rescale = os.path.join(self.test_outdir, 'output_multiband_rescale0-5.tif')
            with rasterio.open(os.path.normpath(out_rescale), 'w', **out_meta) as out:
                out.write_band(1, rescaled0_5)

        self.assertEqual(0, int(np.nanmin(rescaled0_255)), )
        self.assertEqual(255, int(np.nanmax(rescaled0_255)),)

        self.assertEqual(0, int(np.nanmin(rescaled0_5)),)
        self.assertEqual(5, int(np.nanmax(rescaled0_5)))

    def test_normalise(self):
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        with rasterio.open(os.path.normpath(self.singletif)) as src:
            norm = raster_ops.normalise(src)
            out_meta = src.meta.copy()

        if KEEP_TEST_OUTPUTS:
            out_norm = os.path.join(self.test_outdir,'output_singletif_normalised.tif')
            with rasterio.open(os.path.normpath(out_norm), 'w', **out_meta) as out:
                out.write_band(1, norm)

        self.assertEqual(3.3085429668426514, float(np.nanmax(norm)), )
        self.assertEqual( -2.758878707885742, float(np.nanmin(norm)),)

    def test_focalStatisticsSingleBand(self):
        with rasterio.open(os.path.normpath(self.singletif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanmean)
            arr = arr.astype(np.float32)
            self.assertEqual(0.0363851934671402, float(np.nanmin(arr)) )
            self.assertEqual( 6.991053581237793, float(np.nanmax(arr)) )
            self.assertEqual('mean5x5_dummy_singleband_94mga54',col, 'Incorrect column name')

            arr, col = raster_ops.focal_statistics(src, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEqual( 2.632201585583971e-06, float(np.nanmin(arr)))
            self.assertEqual(0.5194834470748901, float(np.nanmax(arr)) )
            self.assertEqual('cv3x3_dummy_singleband_94mga54', col, 'Incorrect column name')

            arr, col = raster_ops.focal_statistics(src, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEqual(5.0, float(np.nanmin(arr)),)
            self.assertEqual(81.0, float(np.nanmax(arr)), )
            self.assertEqual('pixelcount9x9_dummy_singleband_94mga54', col, 'Incorrect column name')

    def test_focalStatisticsMultiBand(self):
        with rasterio.open(os.path.normpath(self.multitif)) as src:
            arr,col = raster_ops.focal_statistics(src, size=5, function=np.nanstd)
            arr = arr.astype(np.float32)
            self.assertEqual( 1.1673586413962767e-05, float(np.nanmin(arr)))
            self.assertEqual(0.3216133415699005, float(np.nanmax(arr)) )
            self.assertEqual('std5x5bd1_dummy_3band_94mga54', col, 'Incorrect column name')

            arr, col = raster_ops.focal_statistics(src, band_num=3, size=3, function=raster_ops.nancv)
            arr = arr.astype(np.float32)
            self.assertEqual(-34.32146072387695, float(np.nanmin(arr)) )
            self.assertEqual(517.6610717773438, float(np.nanmax(arr)) )
            self.assertEqual('cv3x3bd3_dummy_3band_94mga54', col, 'Incorrect column name')

            arr, col = raster_ops.focal_statistics(src, band_num=2, size=9, function=raster_ops.pixelcount)
            arr = arr.astype(np.float32)
            self.assertEqual(5.0, float(np.nanmin(arr)))
            self.assertEqual(81.0, float(np.nanmax(arr)))
            self.assertEqual('pixelcount9x9bd2_dummy_3band_94mga54', col, 'Incorrect column name')
