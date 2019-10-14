import unittest
import utils
import os
import sys
import re
import modeller
import subprocess
from modeller import soap_protein_od
from modeller.terms import energy_term

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.soap

# Mock the soap_protein_od.Scorer class, so a) we can test it without having
# the SOAP potentials installed (even if we do have them, they are slow to
# read in and use a lot of memory) and b) we can make sure any exceptions
# are handled properly by the soap script
class MockScorer(energy_term):
    name = 'MockScorer'
    def __init__(self):
        self.count = 0
    def _assess(self, atmsel, schedule_scale=None, **vars):
        self.count += 1
        if self.count == 1:
            return 10., ()
        elif self.count == 2:
            return 20., ()
        elif self.count == 3:
            raise modeller.ModellerError("SOAP error")

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to soap"""
        for args in (['x'],):
            out = utils.check_output(['cryptosite', 'soap'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output([sys.executable, '-m',
                                     'cryptosite.soap'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_soap_module(self):
        """Test soap script as a module"""
        with utils.temporary_working_directory() as tmpdir:
            for pdb in ('test.pdb', 'pm.pdb0.pdb', 'pm.pdb1.pdb',
                        'pm.pdb2.pdb', 'pm.pdb.B10010002.pdb'):
                fname = os.path.join(tmpdir, pdb)
                with open(fname, 'w') as fh:
                    fh.write("""
ATOM      1  N   CYS A   1      18.511  -1.416  15.632  1.00  6.84           C
""")
            with utils.mocked_object(soap_protein_od, 'Scorer', MockScorer):
                cryptosite.soap.soap_score()
            with open('SnapList.txt') as fh:
                data = fh.readlines()
            # don't know glob order
            self.assertEqual(len(data), 2)
            self.assertTrue(re.match('pm\.pdb\d\.pdb\t10\.0', data[0]))
            self.assertTrue(re.match('pm\.pdb\d\.pdb\t20\.0', data[1]))

    def test_soap_subproc(self):
        """Test soap script as a subprocess"""
        env = os.environ.copy()
        pypath = os.path.join(TOPDIR, 'test', 'mock')
        env['PYTHONPATH'] = pypath + ':' + env.get('PYTHONPATH', '')
        subprocess.check_call(['cryptosite', 'soap'], env=env)

if __name__ == '__main__':
    unittest.main()
