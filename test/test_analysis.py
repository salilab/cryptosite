import unittest
import utils
import os
import sys
import re
import shutil
import subprocess
import modeller

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.analysis

def touch(fname):
    with open(fname, 'w') as fh:
        pass

class MockAtom(object):
    def __init__(self, name, x, y, z):
        self.name = name
        self.x, self.y, self.z = x, y, z

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to analysis"""
        for args in ([],):
            out = utils.check_output(['cryptosite', 'analysis'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output([sys.executable, '-m',
                                     'cryptosite.analysis'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_get_coordinates_sc(self):
        """Test get_coordinates_sc() function"""
        e = modeller.environ()
        m = modeller.model(e)
        coord = cryptosite.analysis.get_coordinates_sc(m,
                os.path.join(TOPDIR, 'test', 'input', 'test_coord.pdb'))
        self.assertEqual(len(coord), 4)
        # First residue is a GLY with no CA -> no coordinates
        self.assertEqual(coord[0], None)
        # Second residue is a GLY with CA -> coordinates are those of
        # the GLY (0,0,0)
        self.assertAlmostEqual(coord[1].x, 0., places=1)
        self.assertAlmostEqual(coord[1].y, 0., places=1)
        self.assertAlmostEqual(coord[1].z, 0., places=1)
        # Third residue is a MET with no sidechain -> no coordinates
        self.assertEqual(coord[2], None)
        # Fourth residue is a MET with a sidechain -> mean coordinates returned
        self.assertAlmostEqual(coord[3].x, 5., places=1)
        self.assertAlmostEqual(coord[3].y, 10., places=1)
        self.assertAlmostEqual(coord[3].z, 15., places=1)

    def test_get_distance(self):
        """Test get_distance() function"""
        a = MockAtom('CA', 0., 0., 0.)
        b = MockAtom('CA', 1., 1., 1.)
        self.assertEqual(cryptosite.analysis.get_distance(a, b), 3.)
        self.assertEqual(cryptosite.analysis.get_distance(a, None), None)
        self.assertEqual(cryptosite.analysis.get_distance(None, b), None)
        self.assertEqual(cryptosite.analysis.get_distance(None, None), None)

    def test_get_qi_bad_num_res(self):
        """Test get_qi with mismatched numbers of residues"""
        class MockModel(object):
            def read(self, file):
                pass
        m = MockModel()
        m.residues = []
        self.assertRaises(ValueError, cryptosite.analysis.get_qi, m, 50,
                          None, None, None)

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

    def test_get_qioft(self):
        """Test get_qioft() function"""
        inputs = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_directory() as tmpdir:
            subdir = os.path.join(tmpdir, 'XXX.pdb_14')
            os.mkdir(subdir)
            for f in ('list', 'pm.pdb.B10010001.pdb', 'pm_XXX.pdb'):
                shutil.copy(os.path.join(inputs, f), subdir)

            cryptosite.analysis.get_qioft(tmpdir)
            with open(os.path.join(subdir, 'qioft_pm_XXX.pdb_11sc.dat')) as fh:
                e = fh.readlines()
            self.assertEqual(e, ['0.2205 0.3422 0.5290 0.3791 0.2468 '
                         '0.4816 0.3841 0.5149 0.4896 0.0971 0.3187 0.0189 \n'])

if __name__ == '__main__':
    unittest.main()
