import glob
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

this_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data')

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

class TestVesperControl(unittest.TestCase):
    def test_updated_keys(self):
        vc = VesperControl(title='This is a test tile')
        self.assertEqual(1, len(vc.updated_keys()))
        vc.update({'datfil': 'MyDataFile.txt',
                   'outdir': 'c:/data/temp'})

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
                   'outdir': "'c:/data/temp'",
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
                           outdir='c:/data/temp',
                           modtyp='Exponential',
                           minpts=150,
                           maxpts=200)

        with self.assertRaises(ValueError) as msg:
            vc.write_to_file(os.path.join('c:/nonexistent/path', 'test_control.txt'))
        self.assertEqual('Output folder c:/nonexistent/path does not exist',
                         str(msg.exception))

        out_file = os.path.join(tempfile.gettempdir(), 'test_control.txt')
        vc.write_to_file(out_file)

        self.assertTrue(out_file)

        with open(out_file) as r_file:
            data = r_file.read()

        self.assertIn("datfil='MyDataFile.txt'", data)
        self.assertIn("outdir='c:/data/temp'", data)
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

    def test1a_CreateControlHighDensity_block_size(self):
        file_csv = os.path.realpath(this_dir + "/area2_high_trimmed.csv")
        grid_file = os.path.realpath(this_dir + "/rasters/area2_5m_blockgrid_v.txt")
        data_col = r'Yld Mass(Dry)(tonne/ha)'

        csv_desc = CsvDescribe(file_csv)
        df_csv = csv_desc.open_pandas_dataframe()

        # check using block_size argument - backwards compatible
        file_bat, file_ctrl = prepare_for_vesper_krige(df_csv, data_col, grid_file, TEMPDIR,
                                                       control_textfile='test_high_5m_control.txt',
                                                       block_size=30, coord_columns=[],
                                                       epsg=28354)
        if os.path.exists(kriging_ops.vesper_exe):
            self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper/Do_Vesper.bat')))
        else:
            self.assertEqual('', file_bat)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_vesperdata.csv')))

        src_df = pd.read_csv(os.path.realpath(this_dir +'/VESPER/high_5m_vesperdata.csv'))
        test_df = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_vesperdata.csv'))

        pd.testing.assert_frame_equal( src_df,test_df)


    def test1_CreateControlHighDensity_VesperControlClass(self):
        # check using VesperControl class
        file_csv = os.path.realpath(this_dir + "/area2_high_trimmed.csv")
        grid_file = os.path.realpath(this_dir + "/rasters/area2_5m_blockgrid_v.txt")
        data_col = r'Yld Mass(Dry)(tonne/ha)'

        csv_desc = CsvDescribe(file_csv)
        df_csv = csv_desc.open_pandas_dataframe()

        vc = VesperControl()
        vc.update(xside=30, yside=30)
        global g_ctrl_file
        file_bat, g_ctrl_file = prepare_for_vesper_krige(df_csv, data_col, grid_file, TEMPDIR,
                                                       control_textfile='test_high_5m_control.txt',
                                                       coord_columns=[],
                                                       epsg=28354,
                                                       control_options=vc)
        if os.path.exists(kriging_ops.vesper_exe):
            self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper/Do_Vesper.bat')))
        else:
            self.assertEqual('', file_bat)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_vesperdata.csv')))

        src_df = pd.read_csv(os.path.realpath(this_dir +'/VESPER/high_5m_vesperdata.csv'))
        test_df = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'test_high_5m_vesperdata.csv'))

        pd.testing.assert_frame_equal( src_df,test_df)

        with open(os.path.realpath(this_dir + '/VESPER/high_5m_control.txt')) as src_file,\
            open(g_ctrl_file) as test_file:
            self.assertEqual(src_file.readlines()[11:], test_file.readlines()[11:])

    def test2_prepareForVesperKrig_LowDensity(self):
        file_csv = os.path.realpath(this_dir + '/area2_lowdensity_points.csv')
        grid_file = os.path.realpath(this_dir + "/rasters/area2_5m_blockgrid_v.txt")
        data_col = 'Hand_Sample'

        csv_desc = CsvDescribe(file_csv)

        ctrl_para = VesperControl()
        ctrl_para.update({'jpntkrg': 1, 'jlockrg': 0, 'minpts': csv_desc.row_count - 2,
                          'maxpts': csv_desc.row_count, 'jcomvar': 0, 'modtyp': 'Spherical',
                          'iwei': 'no_pairs/variance', 'CO': 92.71, 'C1': 277.9, 'A1': 116.0 })

        file_bat, file_ctrl = prepare_for_vesper_krige(csv_desc.open_pandas_dataframe(),
                                                       data_col, grid_file , TEMPDIR,
                                                       control_textfile='test_low_control.txt',
                                                       control_options=ctrl_para,
                                                       coord_columns=[], epsg=28354)
        if os.path.exists(kriging_ops.vesper_exe):
            self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper/Do_Vesper.bat')))
        else:
            self.assertEqual('', file_bat)

        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_low_control.txt')))
        self.assertTrue(os.path.exists(os.path.join(TEMPDIR, 'Vesper', 'test_low_vesperdata.csv')))

        df_csv = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'test_low_vesperdata.csv'))
        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())

        with open(file_ctrl) as r_file:
            data = r_file.read()

        self.assertIn("outfil='test_low_kriged.txt'", data)
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

        src_df = pd.read_csv(os.path.realpath(this_dir +'/VESPER/low_vesperdata.csv'))
        test_df = pd.read_csv(os.path.join(TEMPDIR, 'Vesper', 'test_low_vesperdata.csv'))

        pd.testing.assert_frame_equal( src_df,test_df)

        with open(os.path.realpath(this_dir + '/VESPER/low_control.txt')) as src_file,\
            open(os.path.join(TEMPDIR, 'Vesper', 'test_low_control.txt')) as test_file:
            self.assertEqual(src_file.readlines()[11:], test_file.readlines()[11:])

    def test3_vesperTextToRaster(self):
        # copy the test data to temp.
        for ea_file in glob.glob(os.path.realpath(this_dir +'/VESPER/high_5m*' )):
            shutil.copy(ea_file, TEMPDIR)

        ctrl_file = os.path.realpath(TEMPDIR +'/high_5m_control.txt')
        # these are requirements so check first
        self.assertTrue(os.path.exists(ctrl_file))

        out_pred_tif, out_se_tif, out_ci_txt = vesper_text_to_raster(ctrl_file, 28354)
        for eaFile in [out_pred_tif, out_se_tif, out_ci_txt]:
            self.assertTrue(os.path.exists(eaFile))

        src_pred_file = os.path.realpath(this_dir + '/VESPER/' +os.path.basename(out_pred_tif))
        src_se_file = os.path.realpath(this_dir + '/VESPER/' + os.path.basename(out_se_tif))
        src_ci_file = os.path.realpath(this_dir + '/VESPER/' + os.path.basename(out_ci_txt))

        with rasterio.open( src_pred_file ) as src_pred,\
             rasterio.open(os.path.normpath(out_pred_tif)) as test_pred:

            self.assertEqual(src_pred.profile,test_pred.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test_pred.crs)

            import numpy as np
            np.testing.assert_array_equal(src_pred.read(), test_pred.read())

        with rasterio.open(src_se_file) as src_pred, \
                rasterio.open(os.path.normpath(out_se_tif)) as test_pred:
            self.assertEqual(src_pred.profile, test_pred.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test_pred.crs)

            import numpy as np
            np.testing.assert_array_equal(src_pred.read(), test_pred.read())

        with open(src_ci_file) as src_file, open(out_ci_txt) as test_file:
            self.assertEqual(src_file.readlines()[:2], test_file.readlines()[:2])


    def test4_RunVesper(self):
        if platform.system() != 'Windows':
            print ('Skipping test4_RunVesper - VESPER only present on Windows')
        else:

            # for this to work this file needs updating first
            # these are requirements so check first
            self.assertTrue(os.path.exists(g_ctrl_file))

            vesper_exe = kriging_ops.vesper_exe

            os.path.join(TEMPDIR, 'Vesper', 'high_kriged.tif')
            try:
                print('Running Vesper, Please wait....')
                run_vesper(g_ctrl_file)
            except IOError as msg:
                self.assertIn('does not exist. Please install and configure for kriging to occur',
                             str(msg))

