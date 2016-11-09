#!/usr/bin/env python

"""Analyze AllosMod results."""

import sys
import os
import glob
import math

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

def _get_all_energies(dirname):
    with open(os.path.join(dirname, 'energy.dat')) as fh:
        for line in fh:
            yield float(line.split()[1])

def get_probability(landscape):
    """Get probability for each model in the given landscape directory.
       get_energy(landscape) needs to have been called first.
       Results are written out into files called p.dat in each directory."""
    import numpy
    # Get stats for entire landscape
    all_energies = []
    for dirname in _get_subdirectories(landscape):
        all_energies.extend(_get_all_energies(dirname))
    all_energies = numpy.array(all_energies)
    emin = numpy.min(all_energies)
    RT = numpy.std(all_energies)
    zpart1 = sum(math.exp(-1. * (e - emin) / RT) for e in all_energies)

    for dirname in _get_subdirectories(landscape):
        with open(os.path.join(dirname, 'p.dat'), 'w') as fh:
            for e in _get_all_energies(dirname):
                p = math.exp(-1. * (e - emin) / RT) / zpart1
                fh.write("%15.13f\n" % p)

def main():
    for landscape in sys.argv[1:]:
        get_energy(landscape)
        get_probability(landscape)

if __name__ == '__main__':
    main()
