import unittest
import utils
import os
import sys
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)
import BMIFeatureParser

class Tests(unittest.TestCase):
    def test_get_cnc(self):
        """Test get_cnc() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        '1abc.pdb')
            res, num = BMIFeatureParser.get_cnc('1abc')
            self.assertEqual(len(res), 8)
            self.assertEqual(num, ('1', '1'))

if __name__ == '__main__':
    unittest.main()
