import unittest
import utils
import os
import sys
import modeller

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
from cryptosite import cleaning

pdb_line = "ATOM      2  CA  ALA A   1     -17.308  -2.623  35.999  " \
           "1.00 25.00           C"

class Tests(unittest.TestCase):

    def test_get_pdb_seq(self):
        """Test get_pdb_seq()"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'test.pdb')
            with open(fname, 'w') as fh:
                fh.write("""
ATOM      1  N   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      2  C   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      3  N   GLY A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      4  C   GLY A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      5  N   TYR B   3      18.511  -1.416  15.632  1.00  6.84           C
ATOM      6  C   TYR B   3      18.511  -1.416  15.632  1.00  6.84           C
""")
            seq = cleaning.get_pdb_seq(fname, 'A')
            self.assertEqual(seq, 'CG')

    def test_detect_invalid_residue_types_ok(self):
        """Test _detect_invalid_residue_types() with OK sequence"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'test.pdb')
            with open(fname, 'w') as fh:
                fh.write(pdb_line + '\n')
            e = modeller.environ()
            m = modeller.model(e, file=fname)
        cleaning._detect_invalid_residue_types(m)

    def test_detect_invalid_residue_types_bad(self):
        """Test _detect_invalid_residue_types() with bad sequence"""
        with utils.temporary_directory() as tmpdir:
            fname = os.path.join(tmpdir, 'test.pdb')
            with open(fname, 'w') as fh:
                fh.write("""
ATOM      1  N   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      2  C   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
ATOM      3  N   HIE A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      4  C   HIE A   2      18.511  -1.416  15.632  1.00  6.84           C
ATOM      5  N   HSD B   3      18.511  -1.416  15.632  1.00  6.84           C
ATOM      6  C   HSD B   3      18.511  -1.416  15.632  1.00  6.84           C
""")
            e = modeller.environ()
            m = modeller.model(e, file=fname)
        self.assertRaises(cleaning.InvalidResiduesError,
                          cleaning._detect_invalid_residue_types, m)

    def test_muscle_align(self):
        """Test muscleAlign()"""
        out = cleaning.muscleAlign('ACGV', 'CG', '1abc', 'A')
        self.assertEqual(out, ('ACGV', '-CG-'))

    def test_get_gaps_single(self):
        """Test get_gaps(), single chain"""
        gaps = cleaning.get_gaps(os.path.join(TOPDIR, 'test', 'input',
                                              'test.ali'))
        self.assertEqual(gaps, ['3:,4:', '10:,12:'])

    def test_get_gaps_multiple(self):
        """Test get_gaps(), multiple chain"""
        gaps = cleaning.get_gaps(os.path.join(TOPDIR, 'test', 'input',
                                              'test-multi.ali'))
        self.assertEqual(gaps, ['3:A,4:A', '10:B,12:B'])

    def test_build_model_no_gaps(self):
        """Test build_model() with no gaps"""
        with utils.temporary_working_directory() as tmpdir:
            with open('alignment.pir', 'w') as fh:
                fh.write(">P1;XXX\n")
                fh.write("structureX:input.pdb:1:A:1:A::::\n")
                fh.write("A*\n")
                fh.write(">P1;xxx_X\n")
                fh.write("sequence:input::::::::\n")
                fh.write("A*\n")
            with open('input.pdb', 'w') as fh:
                fh.write(pdb_line + '\n')
            cleaning.build_model('XXX', ['A'])
            os.unlink('XXX_mdl.pdb')

    def test_build_model_gaps(self):
        """Test build_model() with gaps"""
        import modeller.automodel
        def mocked_loopmodel_make(self):
            self.loop.outputs.append({'failure':None, 'name':'input.pdb',
                                      'Normalized DOPE score':-1.})
        with utils.temporary_working_directory() as tmpdir:
            with open('alignment.pir', 'w') as fh:
                fh.write(">P1;XXX\n")
                fh.write("structureX:input.pdb:1:A:2:A::::\n")
                fh.write("A--A*\n")
                fh.write(">P1;xxx_X\n")
                fh.write("sequence:input::::::::\n")
                fh.write("AAAA*\n")
            with open('input.pdb', 'w') as fh:
                fh.write(pdb_line + '\n')
                fh.write(pdb_line[:25] + '2' + pdb_line[26:] + '\n')
            with utils.mocked_object(modeller.automodel.loopmodel, 'make',
                                     mocked_loopmodel_make):
                cleaning.build_model('XXX', ['A'])
            os.unlink('XXX_mdl.pdb')

if __name__ == '__main__':
    unittest.main()
