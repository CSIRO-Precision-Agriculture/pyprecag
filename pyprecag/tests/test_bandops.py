import os
import shutil
import tempfile
import time
import unittest

from pyprecag.bandops import CalculateIndices, BandMapping

pyFile = os.path.basename(__file__)
this_dir = os.path.abspath(os.path.dirname(__file__))

TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])


class TestBandOps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestBandOps, cls).setUpClass()
        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)

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

    def test_band_mapping(self):
        with self.assertRaises(AttributeError) as msg:
            bm = BandMapping(red=4, green=3, blud=2)
        self.assertEqual("BandMapping has no attribute blud. Allowed keys are "
                "blue, mask, infrared, rededge, green, red", str(msg.exception))

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

        file_image = os.path.realpath(this_dir + "/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif")

        with self.assertRaises(KeyError) as msg:
            ci.calculate('NDVIa', file_image, src_nodata=0)
        self.assertEqual("'NDVIA'", str(msg.exception))
