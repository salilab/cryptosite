import unittest
import utils
import os
import sys
import re
import subprocess
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.mainer_short

# Simply copy pre-built outputs at each step, rather than running BLAST,
# fpocket, etc. (these are tested elsewhere)
def mock_muscle_align(query, sbjct, pdb, chain):
    return ('YVTEPCI', 'YVTEPCI')

def mock_run_blast(pdbchain):
    pass

def mock_parse_blast(infile, pdbchain, seq):
    shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdlA.sqc'),
                'XXXA.sqc')

def mock_build_model(pdb, chains):
    shutil.copy('XXX.pdb', 'XXX_mdl.pdb')

def mock_HydChrSSE(pdb, chain):
    shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdlA.hcs'), '.')

def mock_patchmap_feature(pdb):
    shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdl.pdb.ptm'), '.')

def mock_gather_features(pdb, chain_order):
    shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdlA.feat'), '.')

def mock_res_parser(fil):
    shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'XXX_mdl.bmiftr'), '.')

class Tests(unittest.TestCase):
    def test_main(self):
        """Test of complete run of mainer"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'test.pdb'),
                        'input.pdb')
            with open('param.txt', 'w') as fh:
                fh.write("dummy\ndummy\nA\n")
            with utils.mocked_objects(
                    [(cryptosite.cleaning, 'muscleAlign', mock_muscle_align),
                     (cryptosite.seq_conservation, 'run_blast', mock_run_blast),
                     (cryptosite.seq_conservation, 'parse_blast',
                      mock_parse_blast),
                     (cryptosite.cleaning, 'build_model', mock_build_model),
                     (cryptosite.hyd_chr_sse, 'HydChrSSE', mock_HydChrSSE),
                     (cryptosite.patch_mapper, 'patchmap_feature',
                      mock_patchmap_feature),
                     (cryptosite.bmi_feature_parser, 'gather_features',
                      mock_gather_features),
                     (cryptosite.res_parser_bmi, 'res_parser',
                      mock_res_parser)]):
                cryptosite.mainer_short.main()
            # Verify that AllosMod inputs were created
            os.unlink('XXX/align.ali')
            os.unlink('XXX/XXX.pdb')
            os.unlink('XXX/input.dat')
            os.unlink('XXX/list')

if __name__ == '__main__':
    unittest.main()
