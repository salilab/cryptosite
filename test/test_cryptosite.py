import unittest
import utils
import os
import sys
import re
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

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

    def test_help_command(self):
        """Check cryptosite command help"""
        out = utils.check_output(['cryptosite', 'help', 'gather'])
        self.assertTrue('Gather all feature information' in out)

    def test_unknown_command(self):
        """Check 'cryptosite' with an unknown command"""
        for args in (['bad-command'], ['help', 'bad-command']):
            out = utils.check_output(['cryptosite'] + args, retcode=1)

if __name__ == '__main__':
    unittest.main()
