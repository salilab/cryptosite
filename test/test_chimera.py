import unittest
import utils
import os
import sys
import re
import shutil
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.chimera

class Tests(unittest.TestCase):

    def test_bad(self):
        """Test wrong arguments to chimera"""
        for args in ([], ['x'] * 4):
            out = utils.check_output(['cryptosite', 'chimera'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output([sys.executable, '-m',
                                     'cryptosite.chimera'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_make_chimera_file(self):
        """Test make_chimera_file() function"""
        cryptosite.chimera.make_chimera_file('url1', 'url2', 'out.chimerax')
        with open('out.chimerax') as fh:
            lines = fh.readlines()
        self.assertEqual(lines[-4], 'open_files("url1", "url2")\n')
        os.unlink('out.chimerax')

if __name__ == '__main__':
    unittest.main()
