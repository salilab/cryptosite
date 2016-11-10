import unittest
import utils
import os
import sys
import shutil
import contextlib

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(TOPDIR, 'lib'))
import cryptosite.seq_conservation

@contextlib.contextmanager
def mock_usearch():
    """Make mock usearch binary"""
    subdir = 'mock_usearch'
    os.mkdir(subdir)
    fname = os.path.join(subdir, 'usearch')
    with open(fname, 'w') as fh:
        fh.write("""#!/usr/bin/env python
import sys
outf = sys.argv[4]
with open(outf, 'w') as fh:
    fh.write('\\t'.join(['S', '0', '292', '*', '*', '*', '*', '*',
                         'AH70_12410', '*\\n']))
    fh.write('\\t'.join(['L', '0', '292', '*', '*', '*', '*', '*',
                         'AH70_12410', '*\\n']))
    fh.write('\\t'.join(['H', '0', '292', '99.7', '+', '0', '0', '292M',
                         'EN70_12566', 'AH70_12410\\n']))
    fh.write('\\t'.join(['S', '1', '292', '*', '*', '*', '*', '*',
                         'EX70_12567', '*\\n']))
    fh.write('\\t'.join(['H', '1', '292', '98.2', '+', '0', '0', '292M',
                         'AH70_12410', 'EX70_12567\\n']))
""")
    os.chmod(fname, 0o755)

    oldpath = os.environ['PATH']
    os.environ['PATH'] = subdir + ':' + os.environ['PATH']
    yield
    os.environ['PATH'] = oldpath

def mock_ucluster(ali, cutoff=0.8):
    return {'G7PIG01':2, 'XXXq':2}

class Tests(unittest.TestCase):
    def test_parse_blast(self):
        """Test parse_blast() function"""
        blast_out = os.path.join(TOPDIR, 'test', 'input', 'XXXA.blast')
        with utils.temporary_working_directory() as tmpdir:
            with utils.mocked_object(cryptosite.seq_conservation, 'ucluster',
                                     mock_ucluster):
                cryptosite.seq_conservation.parse_blast(blast_out, 'XXX',
                                                        'AMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL')
            os.unlink('XXX.ali')
            os.unlink('XXX.sqc')

    def test_ucluster(self):
        """Test ucluster() function"""
        with utils.temporary_working_directory() as tmpdir:
            with mock_usearch():
                clusters = cryptosite.seq_conservation.ucluster('dummy.ali')
            self.assertEqual(len(clusters), 4)

if __name__ == '__main__':
    unittest.main()
