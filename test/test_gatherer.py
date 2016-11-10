import unittest
import utils
import os
import sys
import re
import shutil
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.gatherer

class Tests(unittest.TestCase):

    def test_bad(self):
        """Test wrong arguments to gatherer"""
        for args in ([], ['x'] * 4):
            out = utils.check_output(['cryptosite', 'gatherer'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output(['python', '-m',
                                     'cryptosite.gatherer'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_process_directory(self):
        """Test process_directory() function"""
        indir = os.path.join(TOPDIR, 'test', 'input', 'gather-inputs')
        inputs = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(inputs, 'XXX_mdl.bmiftr'), '.')
            cryptosite.gatherer.process_directory(indir)
            os.unlink('XXX.features')
            os.unlink('XXX.am')

if __name__ == '__main__':
    unittest.main()
