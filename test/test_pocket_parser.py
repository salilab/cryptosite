import unittest
import utils
import os
import sys
import re
import subprocess
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.pocket_parser

class Tests(unittest.TestCase):

    def test_get_cnc(self):
        """Test get_cnc() function"""
        res = cryptosite.pocket_parser.get_cnc(os.path.join(TOPDIR, 'test',
                                                      'input', 'test.pdb'),
                                               None)
        self.assertEqual(len(res), 8)
        self.assertEqual(res[('ILE', 9, 'A')], 0.0)

    def test_main(self):
        """Test simple complete run of pocket_parser"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input',
                                     'pm.pdb.B10010001.pdb'), '.')
            with open('SnapList.txt', 'w') as fh:
                fh.write("pm.pdb.B10010001.pdb -100.0\n")
                fh.write("high-energy.pdb -10.0\n")
            subprocess.check_call(['cryptosite', 'pocket_parser'])
            with open('pockets.out') as fh:
                lines = fh.readlines()
        self.assertEqual(len(lines), 13)
        self.assertEqual(lines[1], 'ALA\t1\tA\t0.0\n')
        self.assertEqual(lines[2], 'MET\t2\tA\t0.092\n')

if __name__ == '__main__':
    unittest.main()
