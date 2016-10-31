import unittest
import utils
import os
import sys
import re
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)
import pocket_parser

class Tests(unittest.TestCase):

    def test_get_cnc(self):
        """Test get_cnc() function"""
        res = pocket_parser.get_cnc(os.path.join(TOPDIR, 'test', 
                                                 'input', 'test.pdb'),
                                    None)
        self.assertEqual(len(res), 8)
        self.assertEqual(res[('ILE', 9, 'A')], 0.0)

if __name__ == '__main__':
    unittest.main()
