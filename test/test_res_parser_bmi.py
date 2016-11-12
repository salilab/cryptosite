import unittest
import utils
import os
import sys
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
from cryptosite import res_parser_bmi

class Tests(unittest.TestCase):
    def test_res_parser(self):
        """Test res_parser() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'pm_XXX.pdb'),
                        'XXX_mdl.pdb')
            for inp in ('XXX_mdlA.feat', 'XXX_mdlA.sqc', 'XXX_mdl.sas',
                        'XXX_mdl.pdb.ptm'):
                shutil.copy(os.path.join(TOPDIR, 'test', 'input', inp), '.')
            res_parser_bmi.res_parser('XXX_mdl')
            with open('XXX_mdl.bmiftr') as fh:
                data = sorted(fh.readlines())
            self.assertEqual(len(data), 3)
            self.assertEqual(data[1][:50],
                    'XXX_mdl\tALA\t1\tA\t2.298\t106.4\t-10.77\t0.0\t'
                    'U\t0.17\t0.0\t')

if __name__ == '__main__':
    unittest.main()
