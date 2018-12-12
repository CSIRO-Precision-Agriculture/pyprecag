import logging
import os
import platform
import shutil
import tempfile
import unittest
import pandas as pd
import rasterio
import time
from pyprecag import convert, processing, config
from pyprecag.describe import predictCoordinateColumnNames
from pyprecag.processing import clean_trim_points
try:
    from pyprecag import kriging_ops
    from pyprecag.kriging_ops import prepare_for_vesper_krige, vesper_text_to_raster, run_vesper
except ImportError as e:
    if "Vesper" not in e.message:
        raise e

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

file = os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv")
poly = os.path.realpath(this_dir + "/data/area2_onebox_94mga54.shp")
data_col = r'Yld Mass(Dry)(tonne/ha)'
fileSubName = os.path.join(TmpDir, os.path.splitext(os.path.basename(file))[0])
block_tif = fileSubName + '_block.tif'


class test_KrigingOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_KrigingOps, cls).setUpClass()
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

    @unittest.skipIf(
        platform.system() != 'Windows',
        'Vesper only present on Windows'
    )
    def test1_prepareForVesperKrig_filesExist(self):

        gdfPoints, gdfPtsCrs = convert.convert_csv_to_points(file, coord_columns_epsg=4326, out_epsg=28354)
        outGDF, outCRS = clean_trim_points(gdfPoints, gdfPtsCrs, data_col, fileSubName + '_trimmed.csv',
                                           poly, thin_dist_m=2.5)

        self.assertTrue(os.path.exists(fileSubName + '_trimmed.csv'))
        self.assertEqual(outGDF.crs, {'init': 'epsg:28354', 'no_defs': True})
        self.assertEqual(len(outGDF), 554)

        processing.block_grid(in_shapefilename=poly,
                              pixel_size=5,
                              out_rasterfilename=block_tif,
                              out_vesperfilename=fileSubName + '_block_v.txt',
                              snap=True,
                              overwrite=True)
        global file_ctrl
        file_bat, file_ctrl = prepare_for_vesper_krige(outGDF, data_col,
                                                       fileSubName + '_block_v.txt', TmpDir,
                                                       control_textfile=os.path.basename(fileSubName) + '_control.txt',
                                                       block_size=30, coord_columns=[], epsg=28354)

        self.assertTrue(os.path.exists(os.path.join(TmpDir, r'Vesper\Do_Vesper.bat')))
        self.assertTrue(
            os.path.exists(os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_control.txt')))
        self.assertTrue(
            os.path.exists(os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_vesperdata.csv')))

        dfCSV = pd.read_csv(os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_vesperdata.csv'))
        x_column, y_column = predictCoordinateColumnNames(dfCSV.columns)
        self.assertEqual(x_column.upper(), 'EASTING')
        self.assertEqual(y_column.upper(), 'NORTHING')

    @unittest.skipIf(
        platform.system() != 'Windows',
        'Vesper only present on Windows'
    )
    def test2_vesperTextToRaster(self):
        file_ctrl = os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_control.txt')

        # these are requirements so check first
        self.assertTrue(os.path.exists(file_ctrl))
        self.assertTrue(os.path.exists(block_tif))

        vesper_exe = kriging_ops.vesper_exe
        self.assertTrue(os.path.exists(vesper_exe))
        os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_kriged.tif')
        if not os.path.exists(os.path.join(TmpDir, r'Vesper', os.path.basename(fileSubName) + '_kriged.tif')):
            print('Running Vesper, Please wait....')
            run_vesper(file_ctrl)

        out_PredTif, out_SETif, out_CITxt = vesper_text_to_raster(file_ctrl, 28354)
        for eaFile in [out_PredTif, out_SETif, out_CITxt]:
            self.assertTrue(os.path.exists(eaFile))

        with rasterio.open(os.path.normpath(out_PredTif)) as dataset:
            self.assertEqual(dataset.count, 1)
            self.assertEqual(dataset.width, 48)
            self.assertEqual(dataset.height, 26)
            self.assertEqual(dataset.nodatavals, (-9999.0,))
            self.assertEqual(dataset.dtypes, ('float32',))
            self.assertEqual(dataset.crs, rasterio.crs.CRS.from_epsg(28354))

        with rasterio.open(os.path.normpath(out_SETif)) as dataset:
            self.assertEqual(dataset.count, 1)
            self.assertEqual(dataset.width, 48)
            self.assertEqual(dataset.height, 26)
            self.assertEqual(dataset.nodatavals, (-9999.0,))
            self.assertEqual(dataset.dtypes, ('float32',))

            self.assertEqual(dataset.crs, rasterio.crs.CRS.from_epsg(28354))
