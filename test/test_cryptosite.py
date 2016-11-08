import unittest
import utils
import os
import sys
import re
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.environ['PATH'] = os.path.join(TOPDIR, 'bin') + ':' + os.environ['PATH']
os.environ['PYTHONPATH'] = os.path.join(TOPDIR, 'lib') + ':' \
                           + os.environ.get('PYTHONPATH', '')
sys.path.append(os.path.join(TOPDIR, 'lib'))

class Tests(unittest.TestCase):

    def test_no_args(self):
        """Check 'cryptosite' with no arguments"""
        out = utils.check_output(['cryptosite'])
        self.assertTrue("Use 'cryptosite help' for help" in out)

    def test_help(self):
        """Check 'cryptosite help'"""
        for args in (['help'], ['help', 'help']):
            out = utils.check_output(['cryptosite'] + args)
            self.assertTrue('Get help on using' in out)

    def test_unknown_command(self):
        """Check 'cryptosite' with an unknown command"""
        for args in (['bad-command'], ['help', 'bad-command']):
            out = utils.check_output(['cryptosite'] + args, retcode=1)

if __name__ == '__main__':
    unittest.main()
