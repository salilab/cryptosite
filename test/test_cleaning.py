import unittest
import utils
import os
import sys

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
from cryptosite import cleaning

class Tests(unittest.TestCase):

    def test_get_pdb_seq(self):
        """Test get_pdb_seq()"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'test.pdb')
            with open(fname, 'w') as fh:
                fh.write("""
ATOM      1  N   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      2  C   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      3  N   GLY A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      4  C   GLY A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      5  N   TYR B   3      18.511  -1.416  15.632  1.00  6.84           C
ATOM      6  C   TYR B   3      18.511  -1.416  15.632  1.00  6.84           C
""")
            seq = cleaning.get_pdb_seq(fname, 'A')
            self.assertEqual(seq, 'CG')

    def test_muscle_align(self):
        """Test muscleAlign()"""
        out = cleaning.muscleAlign('ACGV', 'CG', '1abc', 'A')
        self.assertEqual(out, ('ACGV', '-CG-'))

    def test_get_gaps_single(self):
        """Test get_gaps(), single chain"""
        gaps = cleaning.get_gaps(os.path.join(TOPDIR, 'test', 'input',
                                              'test.ali'))
        self.assertEqual(gaps, ['3:,4:', '10:,12:'])

    def test_get_gaps_multiple(self):
        """Test get_gaps(), multiple chain"""
        gaps = cleaning.get_gaps(os.path.join(TOPDIR, 'test', 'input',
                                              'test-multi.ali'))
        self.assertEqual(gaps, ['3:A,4:A', '10:B,12:B'])

if __name__ == '__main__':
    unittest.main()
