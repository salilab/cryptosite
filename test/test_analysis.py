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
            subdir = os.path.join(tmpdir, 'XXX.pdb_14')
            os.mkdir(subdir)
            touch(os.path.join(subdir, 'pm.pdb.B10010001.pdb'))
            touch(os.path.join(subdir, 'pm.pdb.B10020001.pdb'))
            touch(os.path.join(subdir, 'pm.pdb.B99990001.pdb'))
            shutil.copy(os.path.join(inputs, 'pm.pdb.D00000001'), subdir)

            e = cryptosite.analysis.get_energy(tmpdir)
            with open(os.path.join(subdir, 'energy.dat')) as fh:
                e = fh.readlines()
            self.assertEqual(e,
                    ['   1000      13478.55371   0.0280   0.0836       '
                     '3374.72437        499.52283\n',
                     '   1000      13580.25098   0.0279   0.0794       '
                     '3332.75293        493.31027\n'])

    def test_get_probability(self):
        """Test get_probability() function"""
        with utils.temporary_directory() as tmpdir:
            subdir = os.path.join(tmpdir, 'XXX.pdb_14')
            os.mkdir(subdir)
            with open(os.path.join(subdir, 'energy.dat'), 'w') as fh:
                fh.write("""1 13291.22266 0.0279 0.0755 3336.66089 493.88873
2 13228.57910 0.0285 0.0719 3450.11133 510.68152
3 13289.01074 0.0281 0.0776 3391.32886 501.98062
4 13319.61426 0.0280 0.0708 3335.29370 493.68634
""")
            cryptosite.analysis.get_probability(tmpdir)
            with open(os.path.join(subdir, 'p.dat')) as fh:
                lines = [float(x) for x in fh.readlines()]
            self.assertEqual(len(lines), 4)
            self.assertAlmostEqual(lines[0], 0.1, places=1)
            self.assertAlmostEqual(lines[1], 0.7, places=1)

if __name__ == '__main__':
    unittest.main()
