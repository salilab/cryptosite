import unittest
import utils
import os
import sys
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
from cryptosite import hyd_chr_sse

class Tests(unittest.TestCase):
    def test_HydChrSSE(self):
        """Test HydChrSSE() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        '1abc.pdb')
            hyd_chr_sse.HydChrSSE('1abc', 'A')
            with open('1ab_mdlA.hcs') as fh:
                data = fh.read()
            self.assertEqual(data,
"""    1 \tA\t-\t0.17\t0
    2 \tY\t-\t-0.94\t0
    3 \tV\t-\t0.07\t0
    4 \tI\t-\t-0.31\t0
    5 \tT\t-\t0.14\t0
    6 \tE\tS\t2.02\t-1
    7 \tP\tT\t0.45\t0
    8 \tC\tT\t-0.24\t0
    9 \tI\t-\t-0.31\t0
""")

if __name__ == '__main__':
    unittest.main()
