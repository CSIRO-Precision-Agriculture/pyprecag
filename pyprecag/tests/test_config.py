import tempfile
import time
import unittest
from pyprecag.config import *

pyFile = os.path.basename(__file__)
this_dir = os.path.abspath(os.path.dirname(__file__))

TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])


class TestConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestConfig, cls).setUpClass()

        global testFailed
        testFailed = False

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_get_config_key(self):
        self.assertEqual(get_config_key('crsLookupURL'),u'http://prj2epsg.org/search.json')

    def test_readconfig(self):
        config = read_config()
        self.assertTrue(isinstance(config, OrderedDict))
        self.assertEqual(config['crsLookupURL'],u'http://prj2epsg.org/search.json')
