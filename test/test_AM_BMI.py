import unittest
import utils
import os
import sys
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)
import AM_BMI

class Tests(unittest.TestCase):
    def test_get_sas(self):
        """Test get_sas() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        '1abc.pdb')
	    out = AM_BMI.get_sas('1abc', 1.4)
            self.assertEqual(len(out), 60)

if __name__ == '__main__':
    unittest.main()
