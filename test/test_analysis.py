import unittest
import utils
import os
import sys
import re
import shutil
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
import cryptosite.analysis

def touch(fname):
    with open(fname, 'w') as fh:
        pass

class Tests(unittest.TestCase):

    def test_get_energy(self):
        """Test get_energy() function"""
        inputs = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_directory() as tmpdir:
            touch(os.path.join(tmpdir, 'pm.pdb.B10010001.pdb'))
            touch(os.path.join(tmpdir, 'pm.pdb.B10020001.pdb'))
            touch(os.path.join(tmpdir, 'pm.pdb.B99990001.pdb'))
            shutil.copy(os.path.join(inputs, 'pm.pdb.D00000001'), tmpdir)

            e = cryptosite.analysis.get_energy(tmpdir)
            self.assertEqual(e,
                    ['   1000      13478.55371   0.0280   0.0836       '
                     '3374.72437        499.52283\n',
                     '   1000      13580.25098   0.0279   0.0794       '
                     '3332.75293        493.31027\n'])

if __name__ == '__main__':
    unittest.main()
