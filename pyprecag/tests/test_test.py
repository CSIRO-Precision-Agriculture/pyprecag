# An empty test to test the unittest framework as a part of the GitHib Action CI/CD
import unittest
import os
import shutil
import tempfile
import time

from pyprecag.tests import setup_folder, KEEP_TEST_OUTPUTS
from pyprecag.bandops import CalculateIndices, BandMapping

PY_FILE = os.path.basename(__file__)
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data', 'rasters')

class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    print("PY_FILE=",PY_FILE)
    print("THIS_DIR=",THIS_DIR)
    unittest.main()
