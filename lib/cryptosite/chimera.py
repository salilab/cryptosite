#!/usr/bin/env python

"""Make a Chimera script to display the output."""

from __future__ import print_function, absolute_import
import optparse
import os
import cryptosite.config

def make_chimera_file(pdb_url, feature_url, chimera_file):
    with open(os.path.join(cryptosite.config.datadir, 'script.chimerax')) as fh:
        chimera_session = fh.read()

    chimera_session += """
open_files("%s", "%s")
]]>
</py_cmd>
</ChimeraPuppet>""" % (pdb_url, feature_url)

    with open(chimera_file, 'w') as out:
        out.write(chimera_session)

def parse_args():
    usage = """%prog [opts] <pdb_url> <feature_url> <chimera_file>

Make a Chimera script to display the output.

<pdb_url> should be the URL for the produced PDB file (XXX.pol.pred.pdb) and
<feature_url> the URL for the produced features file (XXX.pol.pred).

A Chimera script is generated as <chimera_file> to display the CryptoSite
output.
"""
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return args

def main():
    pdb_url, feature_url, chimera_file = parse_args()
    make_chimera_file(pdb_url, feature_url, chimera_file)

if __name__ == '__main__':
    main()
