import os
import shutil
import tempfile
import time
import unittest
from pyprecag.tests import setup_folder, KEEP_TEST_OUTPUTS
from pyprecag.bandops import CalculateIndices, BandMapping

PY_FILE = os.path.basename(__file__)
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data', 'rasters')

TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])


class TestBandOps(unittest.TestCase):
    failedTests = []
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestBandOps, cls).setUpClass()
        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print ('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_band_mapping(self):
        with self.assertRaises(AttributeError) as msg:
            bm = BandMapping(red=4, green=3, blud=2)
        self.assertEqual("BandMapping has no attribute blud. Allowed keys are "
                "blue, green, infrared, mask, red, rededge", str(msg.exception))

        self.assertDictEqual(BandMapping(),
                             {'blue': 0, 'mask': 0, 'infrared': 0, 'rededge': 0, 'green': 0,
                              'red': 0})

        mydict = {'rededge': 1, 'red': 3, 'infrared': 4, 'mask': 5, }
        self.assertDictEqual(BandMapping(**mydict),
                             {'blue': 0, 'mask': 5, 'infrared': 4, 'rededge': 1, 'green': 0,
                              'red': 3})

        bm = BandMapping(red=3, green=2)
        self.assertEqual(bm.allocated_bands(), [2, 3])

        bm.update(infrared=4, mask=5)
        self.assertDictEqual(bm, {'blue': 0, 'mask': 5, 'infrared': 4, 'rededge': 0, 'green': 2,
                                  'red': 3})

        del bm['mask']
        self.assertDictEqual(bm, {'blue': 0, 'mask': 0, 'infrared': 4, 'rededge': 0, 'green': 2,
                                  'red': 3})

    def test_calculateIndices(self):

        self.assertEqual(CalculateIndices(red=3, infrared=4, mask=5).valid_indices(),
                         ['NDVI', 'PCD'])
        ci = CalculateIndices()
        self.assertDictEqual(ci.band_map,
                             {'blue': 0, 'mask': 0, 'infrared': 0, 'rededge': 0, 'green': 0,
                              'red': 0})

        ci.band_map.update({'rededge': 1, 'green': 2, 'mask': 5, 'red': 3, 'infrared': 4})
        self.assertDictEqual(ci.band_map,
                             {'blue': 0, 'green': 2, 'infrared': 4, 'mask': 5, 'red': 3,
                              'rededge': 1})
        self.assertEqual(ci.band_map.allocated_bands(), [1, 2, 3, 4, 5])
        self.assertEqual(ci.valid_indices(), ['NDVI', 'PCD', 'GNDVI', 'NDRE', 'CHLRE'])

        bm = BandMapping(red=3, green=2, infrared=4)
        ci.band_map = bm
        self.assertEqual(ci.valid_indices(), ['NDVI', 'PCD', 'GNDVI'])

        file_image = os.path.realpath(os.path.join(THIS_DIR , "area1_rgbi_jan_50cm_84sutm54.tif"))

        with self.assertRaises(KeyError) as msg:
            ci.calculate('NDVIa', file_image, src_nodata=0)
        self.assertEqual("'NDVIA'", str(msg.exception))
