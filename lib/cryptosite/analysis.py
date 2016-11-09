#!/usr/bin/env python

"""Analyze AllosMod results."""

import sys
import os
import glob

def _get_subdirectories(dirname):
    """Get all subdirectories of the given directory, as full paths."""
    return [x for x in glob.glob(os.path.join(dirname, '*'))
            if os.path.isdir(x)]

def get_energy(landscape):
    """Get energy for each model in the given landscape directory.
       This is just the last entry in the Modeller tracefile for each model,
       and is written out into files called energy.dat in each directory."""
    for dirname in _get_subdirectories(landscape):
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
        with open(os.path.join(dirname, 'energy.dat'), 'w') as energy:
            energy.write(''.join(energy_lines[-nmodels:]))

def main():
    for landscape in sys.argv[1:]:
        get_energy(landscape)

if __name__ == '__main__':
    main()
