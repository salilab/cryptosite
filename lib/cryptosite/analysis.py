#!/usr/bin/env python

"""Analyze AllosMod results."""

import sys
import os
import glob

def get_energy(dirname):
    """Get energy for each model in the given directory.
       This is just the last entry in the Modeller tracefile for each model."""
    # Don't count pm.pdb.B99990001.pdb
    nmodels = len(glob.glob(os.path.join(dirname, 'pm.pdb.B[1-8]*.pdb')))
    last_energy_line = None
    energy_lines = []
    with open(os.path.join(dirname, 'pm.pdb.D00000001')) as fh:
        for line in fh:
            if not line.startswith('#'):
                last_energy_line = line
            elif line.startswith('#  Step'):
                energy_lines.append(last_energy_line)
    energy_lines.append(last_energy_line)
    return energy_lines[-nmodels:]

def main():
    for dirname in sys.argv[1:]:
        with open(os.path.join(dirname, 'energy.dat'), 'w') as fh:
            fh.write(''.join(get_energy(dirname)))

if __name__ == '__main__':
    main()
