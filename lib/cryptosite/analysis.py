#!/usr/bin/env python

"""Analyze AllosMod results."""

from __future__ import print_function, absolute_import
import sys
import os
import glob
import math
import optparse

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

class _AveragedAtom(object):
    """A fake atom whose coordinates are the average of other atoms."""
    def __init__(self, atoms):
        self.x = self.y = self.z = 0.
        for a in atoms:
            self.x += a.x
            self.y += a.y
            self.z += a.z
        self.x /= len(atoms)
        self.y /= len(atoms)
        self.z /= len(atoms)

def get_coordinates_sc(m, fname):
    """Get the coordinates of the sidechain of each residue in the file"""
    m.read(file=fname)
    def get_coord(res):
        if res.name == 'GLY':
            if 'CA' in res.atoms:
                return res.atoms['CA']
        else:
            avg_atoms = [a for a in res.atoms
                         if a.name not in ('CA', 'O', 'N', 'C', 'OT')
                         and not a.name.startswith('H')]
            if avg_atoms:
                return _AveragedAtom(avg_atoms)
    return [get_coord(r) for r in m.residues]

def get_distance(ci, cj):
    if ci is not None and cj is not None:
        dx = ci.x - cj.x
        dy = ci.y - cj.y
        dz = ci.z - cj.z
        return dx*dx + dy*dy + dz*dz

def get_distances(coord, rcut):
    rcut2 = rcut * rcut
    dist = {}
    for i in range(len(coord) - 1):
        for j in range(i + 1, len(coord)):
            d = get_distance(coord[i], coord[j])
            if d is not None and d < rcut2:
                dist[(i,j)] = dist[(j,i)] = math.sqrt(d)
    return dist

def get_qi(m, len_coord, dist, template, fh):
    coord = get_coordinates_sc(m, template)
    if len(coord) != len_coord:
        raise ValueError("different numbers of residues")
    qi_cut = [0.] * len(coord)
    nqi_cut = [0] * len(coord)
    for i in range(len(coord)):
        for j in range(len(coord)):
            if abs(i-j) < 2 or (i,j) not in dist:
                continue
            d = get_distance(coord[i], coord[j])
            if d is not None:
                delta = (dist[(i,j)] - math.sqrt(d)) / (abs(j-i) ** 0.15)
                qi_cut[i] += math.exp(-delta * delta * 0.5)
                nqi_cut[i] += 1
    for i in range(len(coord)):
        if nqi_cut[i] > 0:
            qi_cut[i] /= nqi_cut[i]
        else:
            qi_cut[i] = 1.1
        fh.write("%.4f " % qi_cut[i])

def get_qioft(landscape, rcut=11.):
    """Calculate Qi for all models in a landscape."""
    import modeller
    modeller.log.none()
    e = modeller.environ()
    e.io.hetatm = False
    m = modeller.model(e)

    for dirname in _get_subdirectories(landscape):
        with open(os.path.join(dirname, 'list')) as fh:
            temp1 = 'pm_' + fh.readline().strip()
            pm = os.path.join(dirname, temp1)
        models = sorted(glob.glob(os.path.join(dirname, 'pm.pdb.B[1-8]*.pdb')))

        coord = get_coordinates_sc(m, pm)
        dist = get_distances(coord, rcut)
        with open(os.path.join(dirname, 'qioft_%s_%dsc.dat'
                                         % (temp1, int(rcut))), 'w') as fh:
            for model in models:
                get_qi(m, len(coord), dist, model, fh)
                fh.write('\n')

def parse_args():
    usage = """%prog [opts] <landscape ...>

Analyze AllosMod results. The analysis is done for each passed landscape
directory, and the generated statistics are written into .dat files in
each directory. Currently energy and Qi statistics are computed.
"""
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) == 0:
        parser.error("incorrect number of arguments")
    return args

def main():
    landscapes = parse_args()
    for landscape in landscapes:
        get_energy(landscape)
        get_qioft(landscape)

if __name__ == '__main__':
    main()
