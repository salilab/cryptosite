import unittest
import utils
import os
import sys
import re
import subprocess
import shutil
import contextlib

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
import cryptosite.patch_mapper

@contextlib.contextmanager
def mock_patch_dock_lig_score():
    """Make mock patch_dock and ligand_score binaries which run quickly"""
    subdir = 'mock_patch_dock'
    os.mkdir(subdir)
    fname = os.path.join(subdir, 'patch_dock.Linux')
    with open(fname, 'w') as fh:
        fh.write("""#!%s
with open('ligands.ids') as fh:
    ligands = [l.rstrip('\\r\\n') for l in fh.readlines()]
for i, lig in enumerate(ligands):
    with open('test.pdb%%d.res' %% i, 'w') as fh:
        fh.write("receptorPdb     (Str)   test.pdb\\n")
        fh.write("ligandPdb     (Str)   %%s\\n" %% lig)
        fh.write(" # | score | pen.  | Area    | as1   | as2   | as12  | ACE     | hydroph | Energy  |cluster| dist. || Ligand Transformation\\n")
        if i == 13:
            fh.write("   1 |   684 | -0.61 |   72.40 |     0 |     0 |     0 |  -63.66 |    0.00 |    0.00 |     0 | 0.00 || -0.21078 -0.07140 0.71339 1.12228 -7.44875 0.86045\\n")
        fh.write("Best Rmsd Result: 100000000.00 rank -1\\n")
        fh.write("Best Rank Result: 100000000.00 rank 100000\\n")
""" % sys.executable)
    os.chmod(fname, 0o755)

    fname = os.path.join(subdir, 'ligand_score_multiple')
    with open(fname, 'w') as fh:
        fh.write("""#!%s
import sys, os
pdb, ligand, trans = sys.argv[1:]
if not os.path.exists(pdb) or not os.path.exists(ligand):
    raise RuntimeError("No PDB/ligand")
with open(trans) as fh:
    tr = fh.read()
with open('mol2_score.res', 'w') as fh:
    if tr:
        fh.write("Score for ISB.pdb trans 0 is -0.41\\n")
""" % sys.executable)
    os.chmod(fname, 0o755)

    oldpath = os.environ['PATH']
    os.environ['PATH'] = subdir + ':' + os.environ['PATH']
    yield
    os.environ['PATH'] = oldpath

class Tests(unittest.TestCase):

    def test_make_ligand_file(self):
        """Test make_ligand_file() function"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'ligands.ids')
            cryptosite.patch_mapper.make_ligand_file(fname)
            with open(fname) as fh:
                lines = fh.readlines()
            self.assertEqual(len(lines), 16)

    def test_get_ligand_mol2(self):
        """Test get_ligand_mol2() function"""
        fname = cryptosite.patch_mapper.get_ligand_mol2('ACM')
        self.assertTrue(os.path.exists(fname))

    def test_read_ligand_data(self):
        """Test read_ligand_data() function"""
        with utils.temporary_working_directory() as tmpdir:
            cryptosite.patch_mapper.make_ligand_file('ligands.ids')
            xyz, ligands = cryptosite.patch_mapper.read_ligand_data()
            self.assertEqual(len(xyz), 16)
            self.assertEqual(len(ligands), 16)

    def test_patchmap_feature(self):
        """Test patchmap_feature() function"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'), '.')
            with mock_patch_dock_lig_score():
                cryptosite.patch_mapper.patchmap_feature('test')
            with open('test.pdb.ptm') as fh:
                contents = fh.readlines()
            self.assertEqual(len(contents), 8)

    def test_modify_patch_dock_params(self):
        """Test _modify_patch_dock_params() function"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'params.txt')
            with open(fname, 'w') as fh:
                fh.write("""
protLib /chem.lib
ligandSeg 5.0 15.0 0.1 1 1 1 5
clusterParams 0.05 2 1.0 2.0
""")
            cryptosite.patch_mapper._modify_patch_dock_params(fname)
            with open(fname) as fh:
                contents = fh.read()
        self.assertEqual(contents, """
protLib /chem.lib
ligandSeg 5.0 15.0 0.5 1 1 0 5
clusterParams 0.1 3 1.0 2.0
""")

if __name__ == '__main__':
    unittest.main()
