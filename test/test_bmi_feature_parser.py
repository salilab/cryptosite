import unittest
import utils
import os
import sys
import shutil
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
from cryptosite import bmi_feature_parser

def mock_check_call(args):
    pass

class Tests(unittest.TestCase):
    def test_get_cnc(self):
        """Test get_cnc() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        '1abc.pdb')
            res, num = bmi_feature_parser.get_cnc('1abc')
            self.assertEqual(len(res), 8)
            self.assertEqual(num, ('1', '1'))

    def test_get_cnc_pockets(self):
        """Test get_cnc() parsing of pocket information"""
        with utils.temporary_working_directory() as tmpdir:
            os.mkdir('test_out')
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test_out.pdb'),
                        'test_out')
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test_info.txt'),
                        'test_out')
            with utils.mocked_object(subprocess, 'check_call', mock_check_call):
                res, num = bmi_feature_parser.get_cnc('test')
            self.assertAlmostEqual(res[('ALA', 2, 'A')], 0.0, places=1)
            self.assertAlmostEqual(res[('ALA', 1, 'A')], 0.8, places=1)

    def test_gather_features(self):
        """Test gather_features() function"""
        def mocked_get_cnc(apo):
            residues = {}
            for i, res in enumerate(('ALA', 'MET', 'GLU', 'ASN', 'PHE', 'GLN',
                                     'LYS', 'VAL', 'GLU', 'LYS', 'ILE', 'GLY')):
                residues[(res, i + 1, 'A')] = 42.
            return residues, ('1', '1')
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'pm_XXX.pdb'),
                        'XXX_mdl.pdb')
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdlA.sqc'),
                        '.')
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdlA.hcs'),
                        '.')
            with utils.mocked_object(bmi_feature_parser, 'get_cnc',
                                     mocked_get_cnc):
                bmi_feature_parser.gather_features('XXX_mdl', ['A'])
            with open('XXX_mdlA.feat') as fh:
                data = fh.readlines()
            self.assertEqual(len(data), 96)
            self.assertEqual(data[0][:116],
                     'ATOM      1  N   ALA A   1     -17.647  -2.542  37.435  '
                     '1.00 25.00           N'
                     '\t12.16\t11.0\t42.0\t18.98\t-\t0.17\t0.0\t2.00')

if __name__ == '__main__':
    unittest.main()
