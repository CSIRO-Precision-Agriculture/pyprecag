import logging
import os
import platform
import shutil
import tempfile
import unittest
import pandas as pd
import rasterio
import time

from pyprecag import convert, processing
from pyprecag.describe import predictCoordinateColumnNames
from pyprecag.processing import clean_trim_points
from pyprecag import kriging_ops
from pyprecag.kriging_ops import prepare_for_vesper_krige, vesper_text_to_raster, run_vesper

PY_FILE = os.path.basename(__file__)
TEMPDIR = tempfile.gettempdir()
TEMPDIR = os.path.join(TEMPDIR, os.path.splitext(PY_FILE)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

file_csv = os.path.realpath(this_dir + "/data/area2_yield_ISO-8859-1.csv")
poly = os.path.realpath(this_dir + "/data/area2_onebox_94mga54.shp")
data_col = r'Yld Mass(Dry)(tonne/ha)'
file_sub_name = os.path.join(TEMPDIR, os.path.splitext(os.path.basename(file_csv))[0])
block_tif = file_sub_name + '_block.tif'


class TestKrigingOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestKrigingOps, cls).setUpClass()

        if os.path.exists(TEMPDIR):
            print 'Folder Exists.. Deleting {}'.format(TEMPDIR)
            shutil.rmtree(TEMPDIR)

        os.mkdir(TEMPDIR)

        global test_failed
        test_failed = False

    @classmethod
    def tearDownClass(cls):
        if not test_failed:
            print ('Tests Passed .. Deleting {}'.format(TEMPDIR))
            shutil.rmtree(TEMPDIR)

    def setUp(self):
        self.start_time = time.time()

    def tearDown(self):
        t = time.time() - self.start_time
        print("%s: %.3f" % (self.id(), t))

    def run(self, result=None):
        global test_failed
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            test_failed = True

    @unittest.skipIf(platform.system() != 'Windows', 'Vesper only present on Windows')
    def test1_prepareForVesperKrig_filesExist(self):

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(file_csv, coord_columns_epsg=4326,
                                                                out_epsg=28354)
        out_gdf, out_crs = clean_trim_points(gdf_points, gdf_pts_crs, data_col,
                                             file_sub_name + '_trimmed.csv',
                                             poly, thin_dist_m=2.5)

        self.assertTrue(os.path.exists(file_sub_name + '_trimmed.csv'))
        self.assertEqual({'init': 'epsg:28354', 'no_defs': True}, out_gdf.crs)
        self.assertEqual(554, len(out_gdf))

        processing.block_grid(in_shapefilename=poly,
                              pixel_size=5,
                              out_rasterfilename=block_tif,
                              out_vesperfilename=file_sub_name + '_block_v.txt',
                              snap=True,
                              overwrite=True)
        global file_ctrl
        file_bat, file_ctrl = prepare_for_vesper_krige(out_gdf, data_col,
                                                       file_sub_name + '_block_v.txt', TEMPDIR,
                                                       control_textfile=os.path.basename(
                                                           file_sub_name) + '_control.txt',
                                                       block_size=30, coord_columns=[], epsg=28354)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper\Do_Vesper.bat')))
        self.assertTrue(os.path.exists(
            os.path.join(TEMPDIR, r'Vesper', os.path.basename(file_sub_name) + '_control.txt')))
        self.assertTrue(os.path.exists(
            os.path.join(TEMPDIR, r'Vesper', os.path.basename(file_sub_name) + '_vesperdata.csv')))

        df_csv = pd.read_csv(
            os.path.join(TEMPDIR, r'Vesper', os.path.basename(file_sub_name) + '_vesperdata.csv'))
        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())


    @unittest.skipIf(platform.system() != 'Windows', 'Vesper only present on Windows')
    def test2_vesperTextToRaster(self):
        file_ctrl = os.path.join(TEMPDIR, r'Vesper',
                                 os.path.basename(file_sub_name) + '_control.txt')

        # these are requirements so check first
        self.assertTrue(os.path.exists(file_ctrl))
        self.assertTrue(os.path.exists(block_tif))

        vesper_exe = kriging_ops.vesper_exe
        self.assertTrue(os.path.exists(vesper_exe))
        os.path.join(TEMPDIR, r'Vesper', os.path.basename(file_sub_name) + '_kriged.tif')
        if not os.path.exists(os.path.join(TEMPDIR, r'Vesper',
                                           os.path.basename(file_sub_name) + '_kriged.tif')):
            print('Running Vesper, Please wait....')
            run_vesper(file_ctrl)

        out_pred_tif, out_se_tif, out_ci_txt = vesper_text_to_raster(file_ctrl, 28354)
        for eaFile in [out_pred_tif, out_se_tif, out_ci_txt]:
            self.assertTrue(os.path.exists(eaFile))

        with rasterio.open(os.path.normpath(out_pred_tif)) as dataset:
            self.assertEqual(1, dataset.count, 1)
            self.assertEqual(47, dataset.width, 47)
            self.assertEqual(25, dataset.height, 25)
            self.assertEqual((-9999.0,), dataset.nodatavals)
            self.assertEqual(('float32',), dataset.dtypes)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs,)

        with rasterio.open(os.path.normpath(out_se_tif)) as dataset:
            self.assertEqual(1, dataset.count, 1)
            self.assertEqual(47, dataset.width, 47)
            self.assertEqual(25, dataset.height, 25)
            self.assertEqual((-9999.0,), dataset.nodatavals)
            self.assertEqual(('float32',), dataset.dtypes)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs, )
