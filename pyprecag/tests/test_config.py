import tempfile
import time
import unittest
from pyprecag.config import *

PY_FILE = os.path.basename(__file__)
TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])


class TestConfig(unittest.TestCase):
    failedTests = []
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestConfig, cls).setUpClass()

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_get_config_key(self):
        self.assertEqual(get_config_key('local_projected_epsg'), [28349, 28350, 28351, 28352, 28353, 28354, 28355, 28356,2193])

    def test_readconfig(self):
        config = read_config()
        self.assertTrue(isinstance(config, OrderedDict))
        self.assertEqual(config['local_projected_epsg'], [28349, 28350, 28351, 28352, 28353, 28354, 28355, 28356,2193])
