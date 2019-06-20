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
from pyprecag.describe import predictCoordinateColumnNames, CsvDescribe
from pyprecag.processing import clean_trim_points
from pyprecag import kriging_ops
from pyprecag.kriging_ops import prepare_for_vesper_krige, vesper_text_to_raster, run_vesper, \
    VesperControl

PY_FILE = os.path.basename(__file__)
TEMPDIR = tempfile.gettempdir()
TEMPDIR = os.path.join(TEMPDIR, os.path.splitext(PY_FILE)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

file_csv = os.path.realpath(this_dir + "/data/area2_yield_ISO-8859-1.csv")
poly = os.path.realpath(this_dir + "/data/area2_onebox_94mga54.shp")
data_col = r'Yld Mass(Dry)(tonne/ha)'


class TestVesperControl(unittest.TestCase):
    def test_updated_keys(self):
        vc = VesperControl(title='This is a test tile')
        self.assertEqual(1, len(vc.updated_keys()))
        vc.update({'datfil': 'MyDataFile.txt',
                   'outdir': 'c:\data\temp'})

        self.assertEqual(3, len(vc.updated_keys()))

    def test_value_errors(self):
        vc = VesperControl(title='This is a test tile')

        # test the ValueErrors -------------------------------------------
        # test invalid option on the value
        with self.assertRaises(ValueError) as msg:
            vc.update(jpntkrg=4)
        self.assertEqual("4 is an invalid option for jpntkrg. Options are {'Punctual': 1, "
                         "'Block': 0, 'Point': 1}", str(msg.exception))

        # test in valid value string
        with self.assertRaises(ValueError) as msg:
            vc.update(jpntkrg='Pnt')
        self.assertEqual("Pnt is an invalid option for jpntkrg. Options are {'Punctual': 1, "
                         "'Block': 0, 'Point': 1}", str(msg.exception))

        # test string instead of number
        with self.assertRaises(ValueError) as msg:
            vc.update(icol_x='1')
        self.assertEqual("value for key icol_x should be a int - Got str", str(msg.exception))

        # make sure lookups work
        vc.update(modtyp='Matern')
        self.assertEqual(7, vc['modtyp'])

        vc.update({'datfil': 'MyDataFile.txt',
                   'outdir': "'c:\\data\\temp'",
                   'modtyp': 'Exponential',
                   'minpts': 150,
                   'maxpts': 200})

    def test_key_errors(self):
        vc = VesperControl(title='This is a test tile')

        # test the ValueErrors -------------------------------------------
        # test invalid option on the value
        with self.assertRaises(AttributeError) as msg:
            vc.update(modtype=4)
        self.assertEqual("VesperControl has no attribute modtype", str(msg.exception))

        vc.update(title='this is a test')
        vc.update(modtyp=4)

    def test_write_to_file(self):
        vc = VesperControl(datfil="MyDataFile.txt",
                           outdir='c:\\data\\temp',
                           modtyp='Exponential',
                           minpts=150,
                           maxpts=200)

        with self.assertRaises(ValueError) as msg:
            vc.write_to_file(os.path.join(r'c:\nonexistent\path', 'test_control.txt'))
        self.assertEqual(r'Output folder c:\nonexistent\path does not exist',
                         str(msg.exception))

        out_file = os.path.join(tempfile.gettempdir(), 'test_control.txt')
        vc.write_to_file(out_file)

        self.assertTrue(out_file)

        with open(out_file) as r_file:
            data = r_file.read()

        self.assertIn("datfil='MyDataFile.txt'", data)
        self.assertIn("outdir='c:\\data\\temp'", data)
        self.assertIn("modtyp=2", data)
        self.assertIn("minpts=150", data)


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
    def test1_prepareForVesperKrig_HighDensity(self):

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(file_csv,
                                                                coord_columns_epsg=4326,
                                                                out_epsg=28354)

        out_gdf, out_crs = clean_trim_points(gdf_points, gdf_pts_crs, data_col,
                                             os.path.join(TEMPDIR, 'high_trimmed.csv'),
                                             poly, thin_dist_m=2.5)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'high_trimmed.csv')))
        self.assertEqual({'init': 'epsg:28354', 'no_defs': True}, out_gdf.crs)
        self.assertEqual(554, len(out_gdf))

        processing.block_grid(in_shapefilename=poly,
                              pixel_size=5,
                              out_rasterfilename=os.path.join(TEMPDIR, 'high_block.tif'),
                              out_vesperfilename=os.path.join(TEMPDIR, 'high_block_v.txt'),
                              snap=True,
                              overwrite=True)

        # check using block_size argument - backwards compatible
        global file_ctrl
        file_bat, file_ctrl = prepare_for_vesper_krige(out_gdf, data_col,
                                                       os.path.join(TEMPDIR, 'high_block_v.txt'),
                                                       TEMPDIR,
                                                       control_textfile='high_control.txt',
                                                       block_size=30, coord_columns=[],
                                                       epsg=28354)
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper\Do_Vesper.bat')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'high_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'high_vesperdata.csv')))

        df_csv = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'high_vesperdata.csv'))

        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())
        self.assertIn('EN_EPSG', df_csv.columns)

        with open(file_ctrl) as r_file:
            data = r_file.read()

        self.assertIn("outfil='high_kriged.txt'", data)
        self.assertIn("outdir=''", data)
        self.assertIn("jpntkrg=0", data)
        self.assertIn("jlockrg=1", data)
        self.assertIn("minpts=90", data)
        self.assertIn("maxpts=100", data)
        self.assertIn("modtyp=2", data)
        self.assertIn("jcomvar=1", data)
        self.assertIn("iwei=1", data)
        self.assertIn("xside=30", data)
        self.assertIn("yside=30", data)

        # check using VesperControl class
        vc = VesperControl()
        vc.update(xside=50, yside=50)
        file_bat, file_ctrl = prepare_for_vesper_krige(out_gdf, data_col,
                                                       os.path.join(TEMPDIR, 'high_block_v.txt'),
                                                       TEMPDIR,
                                                       control_textfile='high_control.txt',
                                                       coord_columns=[],
                                                       epsg=28354,
                                                       control_options=vc)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper\Do_Vesper.bat')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'high_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'high_vesperdata.csv')))

        df_csv = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'high_vesperdata.csv'))

        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())
        self.assertIn('EN_EPSG', df_csv.columns)

        with open(file_ctrl) as r_file:
            data = r_file.read()

        self.assertIn("outfil='high_kriged.txt'", data)
        self.assertIn("outdir=''", data)
        self.assertIn("jpntkrg=0", data)
        self.assertIn("jlockrg=1", data)
        self.assertIn("minpts=90", data)
        self.assertIn("maxpts=100", data)
        self.assertIn("modtyp=2", data)
        self.assertIn("jcomvar=1", data)
        self.assertIn("iwei=1",data)
        self.assertIn("xside=50", data)
        self.assertIn("yside=50", data)


    @unittest.skipIf(platform.system() != 'Windows', 'Vesper only present on Windows')
    def test2_vesperTextToRaster(self):
        file_ctrl = os.path.join(TEMPDIR, 'Vesper', 'high_control.txt')
        block_tif = os.path.join(TEMPDIR, 'high_block.tif')

        # these are requirements so check first
        self.assertTrue(os.path.exists(file_ctrl))
        self.assertTrue(os.path.exists(block_tif))

        vesper_exe = kriging_ops.vesper_exe
        self.assertTrue(os.path.exists(vesper_exe))
        os.path.join(TEMPDIR, r'Vesper', 'high_kriged.tif')
        if not os.path.exists(os.path.join(TEMPDIR, r'Vesper', 'high_kriged.tif')):
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
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs, )

        with rasterio.open(os.path.normpath(out_se_tif)) as dataset:
            self.assertEqual(1, dataset.count, 1)
            self.assertEqual(47, dataset.width, 47)
            self.assertEqual(25, dataset.height, 25)
            self.assertEqual((-9999.0,), dataset.nodatavals)
            self.assertEqual(('float32',), dataset.dtypes)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs, )

    def test3_prepareForVesperKrig_LowDensity(self):
        file_csv = os.path.realpath(this_dir + '/data/area2_lowdensity_points.csv')

        processing.block_grid(in_shapefilename=poly,
                              pixel_size=5,
                              out_rasterfilename=os.path.join(TEMPDIR, 'low_block.tif'),
                              out_vesperfilename=os.path.join(TEMPDIR, 'low_block_v.txt'),
                              snap=True,
                              overwrite=True)

        csv_desc = CsvDescribe(file_csv)

        ctrl_para = VesperControl()
        ctrl_para.update({'jpntkrg': 1,
                          'jlockrg': 0,
                          'minpts': csv_desc.row_count - 2,
                          'maxpts': csv_desc.row_count,
                          'jcomvar': 0,
                          'modtyp': 'Spherical',
                          'iwei': 'no_pairs/variance',
                          'CO': 92.71,
                          'C1': 277.9,
                          'A1': 116.0,
                          })

        file_bat, file_ctrl = prepare_for_vesper_krige(csv_desc.open_pandas_dataframe(),
                                                       'Hand_Sample',
                                                       os.path.join(TEMPDIR, 'low_block_v.txt')
                                                       , TEMPDIR,
                                                       control_textfile='low_control.txt',
                                                       control_options=ctrl_para,
                                                       coord_columns=[], epsg=28354)
        print('Running Vesper, Please wait....')
        run_vesper(file_ctrl)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper', 'Do_Vesper.bat')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper', 'low_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, r'Vesper', 'low_vesperdata.csv')))

        df_csv = pd.read_csv(os.path.join(TEMPDIR, r'Vesper', 'low_vesperdata.csv'))
        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())

        with open(file_ctrl) as r_file:
            data = r_file.read()

        self.assertIn("outfil='low_kriged.txt'", data)
        self.assertIn("outdir=''", data)
        self.assertIn("jpntkrg=1", data)
        self.assertIn("jlockrg=0", data)
        self.assertIn("minpts=198", data)
        self.assertIn("maxpts=200", data)
        self.assertIn("modtyp=1", data)
        self.assertIn("jcomvar=0", data)
        self.assertIn("iwei=3",data)
        self.assertIn("CO=92.71", data)
        del data

        out_pred_tif, out_se_tif, out_ci_txt = vesper_text_to_raster(file_ctrl, 28354)
        for eaFile in [out_pred_tif, out_se_tif, out_ci_txt]:
            self.assertTrue(os.path.exists(eaFile))

        with rasterio.open(os.path.normpath(out_pred_tif)) as dataset:
            self.assertEqual(1, dataset.count, 1)
            self.assertEqual(47, dataset.width, 47)
            self.assertEqual(25, dataset.height, 25)
            self.assertEqual((-9999.0,), dataset.nodatavals)
            self.assertEqual(('float32',), dataset.dtypes)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs, )

        with rasterio.open(os.path.normpath(out_se_tif)) as dataset:
            self.assertEqual(1, dataset.count, 1)
            self.assertEqual(47, dataset.width, 47)
            self.assertEqual(25, dataset.height, 25)
            self.assertEqual((-9999.0,), dataset.nodatavals)
            self.assertEqual(('float32',), dataset.dtypes)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), dataset.crs, )
