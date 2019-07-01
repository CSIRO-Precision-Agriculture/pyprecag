import os
import tempfile
from unittest import TestCase
import pandas as pd

from pyprecag.table_ops import calculate_strip_stats

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])
this_dir = os.path.abspath(os.path.dirname(__file__))


class TestTableOps(TestCase):
    def test_calculate_strip_stats(self):
        df = pd.DataFrame({'TrialPtID': ["'0-1'", "'0-2'", "'0-3'", "'0-4'", "'0-5'", "'0-6'",
                                         "'0-7'", "'0-8'", "'0-9'"],
                           'N Strip Control': [4.9786, 5.0882, 5.1409, 5.1353, 5.1179, 5.0467,
                                               5.0243, 5.0915, 5.1767],
                           'S Strip Control': [4.6876, 4.7000, 4.6618, 4.7158, 4.7514, 4.8022,
                                               4.8226, 4.8955, 4.9709],
                           'Strip Value': [4.6034, 4.9931, 5.2358, 5.3660, 5.3959, 5.3928,
                                           5.3153, 5.3071, 5.4854]})

        result, col = calculate_strip_stats(df, 'Strip Value', size=5,
                                       control_columns=['N Strip Control', 'S Strip Control'])

        self.assertEqual(9, len(result.columns))
        self.assertEqual('N-S_mean',col)








