import unittest
import utils
import os
import sys
import shutil
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.am_bmi

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to am_bmi"""
        for args in (['x'],):
            out = utils.check_output(['cryptosite', 'am_bmi'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output([sys.executable, '-m',
                                     'cryptosite.am_bmi'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_get_sas(self):
        """Test get_sas() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        '1abc.pdb')
            out = cryptosite.am_bmi.get_sas('1abc', 1.4)
            self.assertEqual(len(out), 60)

    def test_main(self):
        """Test simple complete run of am_bmi"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input',
                                     'pm.pdb.B10010001.pdb'), '.')
            with open('SnapList.txt', 'w') as fh:
                fh.write("pm.pdb.B10010001.pdb -100.0\n")
                fh.write("high-energy.pdb -10.0\n")
            subprocess.check_call(['cryptosite', 'am_bmi'])
            with open('am_features.out') as fh:
                lines = sorted(fh.readlines())
        self.assertEqual(len(lines), 12)
        self.assertEqual(lines[0], 'ALA\t1\tA\t17.328\t12.02\t32.6\t48.0\n')

if __name__ == '__main__':
    unittest.main()
