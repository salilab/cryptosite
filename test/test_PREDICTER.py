import unittest
import utils
import os
import sys
import re
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
sys.path.append(TOPDIR)
import PREDICTER

class Tests(unittest.TestCase):

    def test_get_matrix(self):
        """Test get_matrix() function"""
        fname = os.path.join(TOPDIR, 'test', 'input', 'test.features')
        m, header, indices = PREDICTER.get_matrix(fname, model='final')
        self.assertEqual(len(m), 2)
        self.assertEqual(header, ['SQC', 'PTM', 'CNC_mean_'])
        self.assertEqual(len(indices), 2)

if __name__ == '__main__':
    unittest.main()
