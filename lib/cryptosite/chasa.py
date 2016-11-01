#!/usr/bin/env python

"""
Title: CHASA.py
Author: Pat Fleming, pat.fleming@jhu.edu
Date: Fall, 2004
Requirements: numpy

Updated: Sun Jan 23 11:35:53 EST 2005

Calculate Conditional Hydrophobic Accessible Surface Area (CHASA), backbone
N and O solvation accessiblity and total solvation free energy for polypeptides.

Also identify unsatisfied hydrogen bond donors (N) and acceptors (O) in backbone.

This script calculates the accessible surface area (ASA) conditional upon prior
solvation of the backbone N and O atoms by placing an oxygen atom proximate to
these backbone atoms at acceptable hydrogen bond distance and orientation and
then including these solvation "waters" in the ASA calculation.

The CHASA method and calculation of total solvation energy are described in:
------------------------------------------------------------------------
Patrick J. Fleming, Nicholas C. Fitzkee, Mihaly Mezei, Rajgopal Srinivasan,
and George D. Rose (2004)

"A novel method reveals that solvent water favors polyproline II over
beta-strand conformation in peptides and unfolded proteins:
Conditional Hydrophobic Accessible Surface Area (CHASA)"

Proteins (2005) 14:111-118
------------------------------------------------------------------------
Input:

A standard PDB file is the input. The current version does not deal well with
hydrogens and they should be stripped. (The code deals only with backbone hydrogens
used in current versions of LINUS).

The N-terminal nitrogen is not solvated and thus the CHASA for some atoms proximate to
this atom in cartesian space will be incorrect.

Output:

The output is in the form of a PDB file with results in the occupancy and B-factor
columns and header and TER records.

Example:

COMPND numintHbd     numSolv    num_nonHbd      total_bb_polar
COMPND     84          147           3              111
.
.
.
ATOM      9  N   THR     2      13.719  19.413  27.573  5.00  0.00
ATOM     10  CA  THR     2      13.088  19.661  26.283  0.00  0.00
ATOM     11  C   THR     2      13.561  18.631  25.300  0.00  0.00
ATOM     12  O   THR     2      14.763  18.432  25.121  3.00  1.57
ATOM     13  CB  THR     2      13.527  20.980  25.667  0.00  7.80
ATOM     14  OG1 THR     2      13.307  22.020  26.627  0.00  5.84
ATOM     15  CG2 THR     2      12.704  21.284  24.409  0.00 16.10
.
.
.
ATOM    149  N   ALA    20       9.346  17.206  29.144  0.00  0.49
ATOM    150  CA  ALA    20       8.985  15.930  29.750  0.00  2.64
ATOM    151  C   ALA    20      10.067  15.607  30.760  0.00  3.52
ATOM    152  O   ALA    20      11.193  16.119  30.686 -1.00  1.77
ATOM    153  CB  ALA    20       8.856  14.815  28.714  0.00  2.39
.
.
.
ATOM    499  O   HOH   120       7.787   1.168  13.283  0.00  0.00
ATOM    500  O   HOH   121       9.368  -0.535   6.166  0.00  0.00
ATOM    501  O   HOH   122      22.918  12.605  15.164  0.00  0.00
ATOM    502  O   HOH   123      16.522  12.543  26.647  0.00  0.00
ATOM    503  O   HOH   124      22.421   8.292  28.255  0.00  0.00
ATOM    504  O   HOH   125      11.087  12.195  10.474  0.00  0.00
TER 1424.592  -18.062

where:

        numintHbd = number of internally hydrogen bonded backbone N and O

        numSolv = number of solvation "waters" (max = 5 per backbone polar group)

        num_nonHbd = number of backbone polar groups not satisfied by hydrogen bonding

        total_bb_polar = number of backbone polar groups in the molecule that should
                         be hydrogen bonded ([2 x Number of residues] -1)

        occupancy column (for N and O only):
            (if > 0.00) = number of solvation "waters" accessible to that atom (max = 5)
            (if = -1.00) = this atom is not hydrogen bond satisfied

        B factor column = CHASA in square angstroms

        TER record:
            (first number) = total CHASA for molecule
            (second number) = solvation free energy for molecule

        HOH atoms = solvation "waters" in hydrogen bonding proximity to backbone polar groups
                    (these are the conditional atoms added prior to ASA calculation)

------------------------------------------------------------------------
Much of the code here is modified from the LINUS suite of programs originally
written by Raj Srinivasan (http://roselab.jhu.edu/dist/).

The code for calculation of ASA by the method of Shrake and Rupley was ported from
C code originally written by Nick Fitzkee.

The script requires that numpy be installed.
------------------------------------------------------------------------
"""

USAGE = """
python chasa.py pdbinfile > [pdboutfile]

For non-conditional (traditional) ASA uncomment the following line
near bottom of script:

#   wat_list = None

"""

import sys, os, time
from math import sqrt, pi, sin, cos, acos, fabs
from string import strip
import gzip
import numpy

RADIANS_TO_DEGREES = 180.0/pi
DEGREES_TO_RADIANS = pi/180.0
StringType = type('')
IntType = type(0)

class Linus:
    """Class to describe a (modified) LINUS object

    Instantiation

        An instance is created by calling the class constructor with
        one argument - name of the pdb file.

    Attributes

        o *protein* - instance of LinusMol - the protein that is being
        simulated

        o *hbdpar* - hydrogen bonding energy parameters.

            i.   *use_hbond* - boolean - if true backbone - backbone
                 hydrogen bonds are enabled

            ii.  *use_sidechain_hbond* - boolean - if true sidechain
                 acceptor to backbone donor hydrogen bonds are
                 enabled

            iii. *hbond_distance* - float - optimal distance between
                 backbone donor and backbone acceptor

            iv.  *hbdond_torsion* - float - minium torsion angle between
                 the backbone acceptor and the backbone donor and it's two
                 antecedent atoms( 'C' of previous residue and 'CA' of the
                 residue

            v.   *hbond_probe* - float - distance over which the hydrogen
                 bond energy for backbone donor to backbone acceptor
                 bonds will be scaled to zero from the maximum value

            vi.  *hbond_score_short* - float - backbone donor to backbone
                 acceptor hydrogen bond energy for residues separations of
                 5 or less

            vii. *hbond_score_long* - float - backbone donor to backbone
                 acceptor hydrogen bond energy for residues separations of
                 6 or less

            viii. *sidechain_hbond_distance* - float - optimal distance between
                  backbone donor and sidechain acceptor

            ix.  *sidechain_hbond_torsion* - float - minium torsion angle
                 between the sidechain acceptor and the backbone donor and
                 it's two antecedent atoms ( 'C' of previous residue and 'CA'
                 of the residue)

            x.   *sidechain_hbond_score* - float - energy of a sidechain
                 to backbone hydrogen bond

            xi.  *hbond_winmin* - integer - minimum separation between
                 two residues that can participate in backbone to backbone
                 hydrogen bonds

            xii. *hbond_winmax* - integer - maximum  separation between
                 two residues that can participate in backbone/sidechain
                 acceptor to backbone donor hydrogen bonds

    """

    def __init__(self, pdbfile):
        """Create a LINUS simulation object

            pdbfile - file containing coordinates of the
                      molecule to be simulated in PDB format


        """

        pdbfile = os.path.expanduser(pdbfile)
        self.protein = LinusProtein(pdbfile)
        self.hbdpar = default_hbdpar.copy()
        self.hblist = self.conlist = self.dscnstrlist = None

    def __getattr__(self, name):
        if self.__dict__.has_key(name):
            return self.__dict__[name]
        elif self.hbdpar.has_key(name):
            return self.hbdpar[name]
        else:
            raise AttributeError, 'Undefined attribute %s' % name

    def __setattr__(self, name, value):
        if name in ('pdbfile', 'protein', 'hbdpar', 'hblist'):
            self.__dict__[name] = value
        elif not hasattr(self, name):
            self.__dict__[name] = value
        else:
            raise AttributeError, 'Attribute %s cannot be set directly' % name

    def set_hbond_parameters(self, **kw):
        """Set the hydrogen bond parameters

        Arguments

            series of keyword=value arguments

        Errors

            will raise AttributeError if invalid parameter name is specified
            and TypeError if invalid value for the parameter is supplied

        """

        for key, value in kw.items():
            if self.hbdpar.has_key(key):
                vtype = type(self.hbdpar[key])
                self.hbdpar[key] = converters[vtype](value)
            else:
                raise KeyError, 'Undefined hbond parameter %s' % key

class LinusProtein:
    """Class to describe a protein in LINUS.

    Instantiation

       Instantiation of this class requires one argument - the
       name of the file containing a description of the molecule
       in PDB format.

    Class Attributes

       o *filename* - string - name of file from which the molecule
       information was read

       o *atoms* - list - list of atoms in the molecule. each item
       in this list is an instance of *Atom3d*

       o *num_atoms* - integer - number of atoms in molecule

       o *num_residues* - integer - number of residues in molecule

       o *residue_names* - list - name of each residue in molecule

       o *res_pdb_number* - list - number of each residue in input PDB file

       o *residue_first_atom_indices* - list - index of the the
       first atom of each residue (this is horrendously named)

       o *phi_atoms* - list - contains a reference to the 'C' atom
       of each residue

       o *psi_atoms* - list - contains a reference to the 'O' atom
       of each residue

       o *omega_atoms* - list - contains a reference to the 'CA' atom
       of each residue

       o *chi_torsions* - list of lists - references to the atoms
       corresponding to the various chi torsions of each residue

    Methods

        o *get_atom* - get a reference to an atom in a molecule given
        the serial number of the residue it belongs to and the atoms'
        name

        o *update_coordinates* - recompute cartesian coordinates for
        a molecule after changes to the internal coordinates

        o *radius_of_gyration* - compute the radius of gyration of
        the molecule

        o *reset_conformation* - restore the conformation of the
        molecule to that specified in the input PDB file from which
        the molecule was first created
    """

    def __init__(self, filename):
        """New instance of LinusMol

        Arguments

            o *filename* - string - name of the PDB format file containing
            a description of the molecule

        Result

            New instace of LinusMol
        """

        self.filename = os.path.expanduser(filename)
        protein_from_pdb(self)
        nr = self.num_residues
        fp = self._orig_fp_distance = []
        sp = self._orig_sp_angle = []
        tp = self._orig_tp_torsion = []
        for atom in self.atoms:
            fp.append(atom.fp_distance)
            sp.append(atom.sp_angle)
            tp.append(atom.tp_torsion)

    def get_atom(self, resid, atomname):
        """Get a reference to an atom in a molecule

        Arguments

            o *resid* - integer - serial number of residue in which
            atom is present

            o *atomname* - string - name of atom

        Result

            Returns a reference to the atom.  If an atom with required
            name is not found in the resiude, *ValueError* is raised.
        """
        i1 = self.residue_first_atom_indices[resid]
        i2 = self.residue_first_atom_indices[resid+1]
        for i in range(i1, i2):
            if self.atoms[i].name == atomname:
                return self.atoms[i]
        else:
            raise ValueError, 'No such atom: %s' % atomname


class Atom3d:
    """Class to describe an atom in a molecule

    Attributes

        o *name* - string - name of an atom in PDB format

        o *radius* - float - hard sphere radius of an atom

        o *x*, *y*, *z* - float - cartesian coordinates of atom

        o *first_parent*, *second_parent*, *third_parent* - instances of
        Atom3d - Z-matrix description of an atom

        o *fp_distance*, *sp_angle*, *tp_torsion* - float - distance,
        angle in degrees and torsion angle in degrees made by atom with its
        first, second and third parents

        o *resnum* - integer - index of residue that atom is a part of

    """
    def __init__(self):
        self.name = ''
        self.radius = 0.0
        self.x = self.y = self.z = 0.0
        self.first_parent = self.second_parent = self.third_parent = None
        self.fp_distance = self.sp_angle = self.tp_torsion = 0.0
        self.fp_distance_tmp = self.sp_angle_tmp = self.tp_torsion_tmp = 0.0
        self.resnum = 0

    def __hash__(self):
        return hash(repr(self.resnum) + '_' + self.name)

    def coords(self):
        """Cartesian coordinates of an atom as a list"""
        return [self.x, self.y, self.z]

    def distance(self, other):
        """Euclidean distance between two atoms"""
        return sqrt((self.x - other.x)**2.0 + (self.y - other.y)**2.0 +
                    (self.z - other.z)**2.0)

    def close(self, other, dist=None):
        """check if two atoms are within a specified distance"""

        if dist is None:
            r1 = 0.90*self.radius
            r2 = 0.90*other.radius
            if r1 < 1.0e-6 or r2 < 1.0e-6: return 0
            dist = r1 + r2
        d2 = dist * dist
        dx = (self.x - other.x)**2.0
        if dx >= d2: return 0
        dy = (self.y - other.y)**2.0
        if dy >= d2: return 0
        dz = (self.z - other.z)**2.0
        if dz >= d2: return 0
        return dx + dy + dz < d2

    def int_to_cart(self):
        """generate cartesian coordinates for an atom from it's
        internal coordinates"""


        a2 = self.first_parent
        if a2 is None:
            self.x = self.y = self.z = 0.0
            return

        a3 = self.second_parent
        if a3 is None:
            self.x = self.fp_distance
            self.y = self.z = 0.0
            return


        ang = self.sp_angle * DEGREES_TO_RADIANS
        sina = sin(ang)
        cosa = -cos(ang)

        a4 = self.third_parent
        if a4 is None:
            self.z = 0.0
            self.x = a2.x + self.fp_distance*cosa
            self.y = self.fp_distance * sina
            return

        tor = self.tp_torsion * DEGREES_TO_RADIANS
        sint = sina * sin(tor)
        cost = sina * cos(tor)

        u1x = a3.x - a4.x
        u1y = a3.y - a4.y
        u1z = a3.z - a4.z
        d = 1.0 / sqrt(u1x**2.0 + u1y**2.0 + u1z**2.0)
        u1x = u1x * d
        u1y = u1y * d
        u1z = u1z * d

        u2x = a2.x - a3.x
        u2y = a2.y - a3.y
        u2z = a2.z - a3.z
        d = 1.0 / sqrt(u2x**2.0 + u2y**2.0 + u2z**2.0)
        u2x = u2x * d
        u2y = u2y * d
        u2z = u2z * d

        cosine = u1x*u2x + u1y*u2y + u1z*u2z

        if abs(cosine) < 1.0:
            sine = 1.0/sqrt(1.0 - cosine*cosine)
        else:
            sine = 1.0/sqrt(cosine*cosine - 1.0)

        u3x = sine * (u1y*u2z - u1z*u2y)
        u3y = sine * (u1z*u2x - u1x*u2z)
        u3z = sine * (u1x*u2y - u1y*u2x)

        u4x = cost * (u3y*u2z - u3z*u2y)
        u4y = cost * (u3z*u2x - u3x*u2z)
        u4z = cost * (u3x*u2y - u3y*u2x)

        dist = self.fp_distance
        self.x = a2.x + dist*(u2x*cosa + u4x + u3x*sint)
        self.y = a2.y + dist*(u2y*cosa + u4y + u3y*sint)
        self.z = a2.z + dist*(u2z*cosa + u4z + u3z*sint)

    def commit_internal_coords(self):
        """Save internal coordinates following acceptance of a conformation
        change"""
        self.fp_distance_tmp = self.fp_distance
        self.sp_angle_tmp = self.sp_angle
        self.tp_torsion_tmp = self.tp_torsion


# RESIDUES - description of residues known to LINUS.  Each residue is
# described as a tuple.  The items in tuple are themselves tuples describing
# each atom in the residue.  The fields of the tuple are:

#  0. name of atom
#  1, 2, 3 - first second and third parent of the atom.  The parents
#            are described as a tuple of two values - offset and name.
#            The legal values for offset are 0 or -1. A value of 0
#            means that the parent for the atom is in the same residue
#            and -1 means previous residue
#  4 - Hard sphere radius of atom
#  5 - Radius of atom for use in contact energy calculations.  A value of
#      0 indicates the atom does not participate in contact energy
#      calculations
#  6 - is 1 if atom is a hydrogen bond donor and 0 otherwise
#  7 - is 1 if atom is a hydrogen bond acceptor and 0 otherwise
#  8 - number of bonds separating atom from its own backbone amide
#      nitrogen
#  9 - number of bonds separating atom from the backbone amide nitrogen
#      of the succeeding residue

ACE = (
    (' CA ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.7, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), (-1, ' C  '), (-1, ' CA '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), (-1, ' C  '), 1.35, 0.0, 0, 1, 3, 2)
)

NH2 = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3)
)

NME = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
)

ALA = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.70, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3)
)

ARG = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.7, 2.0, 0, 0, 4, 5),
    (' NE ', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' CB '), 1.35, 0.0, 1, 0, 5, 6),
    (' CZ ', ( 0, ' NE '), ( 0, ' CD '), ( 0, ' CG '), 1.50, 0.0, 0, 0, 6, 7),
    (' NH1', ( 0, ' CZ '), ( 0, ' NE '), ( 0, ' CD '), 1.35, 0.0, 1, 0, 7, 8),
    (' NH2', ( 0, ' CZ '), ( 0, ' NE '), ( 0, ' NH1'), 1.35, 0.0, 1, 0, 7, 8)
)

ASN = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 0.0, 0, 0, 3, 4),
    (' OD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.35, 0.0, 1, 1, 4, 5),
    (' ND2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' OD1'), 1.35, 0.0, 1, 1, 4, 5)
)

ASP = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 0.0, 0, 0, 3, 4),
    (' OD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.35, 0.0, 1, 1, 4, 5),
    (' OD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' OD1'), 1.35, 0.0, 1, 1, 4, 5)
)

CYS = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' SG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 0.0, 1, 1, 3, 4),
)

GLN = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.50, 0.0, 0, 0, 4, 5),
    (' OE1', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' CB '), 1.35, 0.0, 1, 1, 5, 6),
    (' NE2', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' OE1'), 1.35, 0.0, 1, 1, 5, 6)
)

GLU = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.50, 0.0, 0, 0, 4, 5),
    (' OE1', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' CB '), 1.35, 0.0, 1, 1, 5, 6),
    (' OE2', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' OE1'), 1.35, 0.0, 1, 1, 5, 6)
)

GLY = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' HB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
)

HIS = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 0.0, 0, 0, 3, 4),
    (' ND1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.35, 0.0, 1, 1, 4, 5),
    (' CD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' ND1'), 1.50, 0.0, 0, 0, 4, 5),
    (' CE1', ( 0, ' ND1'), ( 0, ' CG '), ( 0, ' CB '), 1.50, 0.0, 0, 0, 5, 6),
    (' NE2', ( 0, ' CD2'), ( 0, ' CG '), ( 0, ' CE1'), 1.35, 0.0, 1, 1, 5, 6),
)

ILE = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG1', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CG2', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' CG1'), 1.7, 2.0, 0, 0, 3, 4),
    (' CD1', ( 0, ' CG1'), ( 0, ' CB '), ( 0, ' CA '), 1.7, 2.0, 0, 0, 4, 5),
)

LEU = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.7, 2.0, 0, 0, 4, 5),
    (' CD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CD1'), 1.7, 2.0, 0, 0, 4, 5),
)

LYS = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.7, 2.0, 0, 0, 4, 5),
    (' CE ', ( 0, ' CD '), ( 0, ' CG '), ( 0, ' CB '), 1.7, 0.0, 0, 0, 5, 6),
    (' NZ ', ( 0, ' CE '), ( 0, ' CD '), ( 0, ' CG '), 1.35, 0.0, 1, 0, 6, 7),
)

MET = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' SD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.80, 0.0, 0, 1, 4, 5),
    (' CE ', ( 0, ' SD '), ( 0, ' CG '), ( 0, ' CB '), 1.7, 2.0, 0, 0, 5, 6),
)

PHE = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.50, 2.0, 0, 0, 3, 4),
    (' CD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.50, 2.0, 0, 0, 4, 5),
    (' CD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CD1'), 1.50, 2.0, 0, 0, 4, 5),
    (' CE1', ( 0, ' CD1'), ( 0, ' CG '), ( 0, ' CB '), 1.50, 2.0, 0, 0, 5, 6),
    (' CE2', ( 0, ' CD2'), ( 0, ' CG '), ( 0, ' CE1'), 1.50, 2.0, 0, 0, 5, 6),
    (' CZ ', ( 0, ' CE1'), ( 0, ' CD1'), ( 0, ' CG '), 1.50, 2.0, 0, 0, 6, 7),
)

PRO = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 0, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 2),
    (' CD ', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.45, 2.0, 0, 0, 1, 1)
)

SER = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' OG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.35, 0.0, 1, 1, 3, 4)
)

THR = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' OG1', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.35, 0.0, 1, 1, 3, 4),
    (' CG2', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' OG1'), 1.7, 2.0, 0, 0, 3, 4)
)

TRP = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.50, 2.0, 0, 0, 3, 4),
    (' CD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.50, 2.0, 0, 0, 4, 5),
    (' CD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CD1'), 1.50, 2.0, 0, 0, 4, 5),
    (' NE1', ( 0, ' CD1'), ( 0, ' CG '), ( 0, ' CB '), 1.35, 0.0, 1, 0, 5, 6),
    (' CE2', ( 0, ' NE1'), ( 0, ' CD1'), ( 0, ' CG '), 1.50, 2.0, 0, 0, 5, 6),
    (' CE3', ( 0, ' CD2'), ( 0, ' CG '), ( 0, ' CE2'), 1.50, 2.0, 0, 0, 5, 6),
    (' CZ2', ( 0, ' CE2'), ( 0, ' NE1'), ( 0, ' CD1'), 1.50, 2.0, 0, 0, 6, 7),
    (' CZ3', ( 0, ' CE3'), ( 0, ' CD2'), ( 0, ' CG '), 1.50, 2.0, 0, 0, 6, 7),
    (' CH2', ( 0, ' CZ3'), ( 0, ' CE3'), ( 0, ' CD2'), 1.50, 2.0, 0, 0, 7, 8),
)

TYR = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG ', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.50, 2.0, 0, 0, 3, 4),
    (' CD1', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CA '), 1.50, 2.0, 0, 0, 4, 5),
    (' CD2', ( 0, ' CG '), ( 0, ' CB '), ( 0, ' CD1'), 1.50, 2.0, 0, 0, 4, 5),
    (' CE1', ( 0, ' CD1'), ( 0, ' CG '), ( 0, ' CB '), 1.50, 2.0, 0, 0, 5, 6),
    (' CE2', ( 0, ' CD2'), ( 0, ' CG '), ( 0, ' CE1'), 1.50, 2.0, 0, 0, 5, 6),
    (' CZ ', ( 0, ' CE1'), ( 0, ' CD1'), ( 0, ' CG '), 1.50, 2.0, 0, 0, 6, 7),
    (' OH ', ( 0, ' CZ '), ( 0, ' CE1'), ( 0, ' CE2'), 1.35, 0.0, 1, 1, 7, 8),
)

VAL = (
    (' N  ', (-1, ' C  '), (-1, ' CA '), (-1, ' O  '), 1.35, 0.0, 1, 0, 0, 3),
    (' CA ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 1.5, 0.0, 0, 0, 1, 2),
    (' C  ', ( 0, ' CA '), ( 0, ' N  '), (-1, ' C  '), 1.50, 0.0, 0, 0, 2, 1),
    (' O  ', ( 0, ' C  '), ( 0, ' CA '), ( 0, ' N  '), 1.34, 0.0, 0, 1, 3, 2),
    (' CB ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.7, 2.0, 0, 0, 2, 3),
    (' H  ', ( 0, ' N  '), (-1, ' C  '), (-1, ' CA '), 0.90, 0.0, 0, 0, 1, 4),
    (' HA ', ( 0, ' CA '), ( 0, ' N  '), ( 0, ' C  '), 1.00, 0.0, 0, 0, 2, 3),
    (' CG1', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' N  '), 1.7, 2.0, 0, 0, 3, 4),
    (' CG2', ( 0, ' CB '), ( 0, ' CA '), ( 0, ' CG1'), 1.7, 2.0, 0, 0, 3, 4),
)
RESIDUES = {
    'ACE': ACE,
    'NH2': NH2,
    'NME': NME,
    'ALA': ALA,
    'ARG': ARG,
    'ASN': ASN,
    'ASP': ASP,
    'CYS': CYS,
    'GLN': GLN,
    'GLU': GLU,
    'GLY': GLY,
    'HIS': HIS,
    'ILE': ILE,
    'LEU': LEU,
    'LYS': LYS,
    'MET': MET,
    'PHE': PHE,
    'PRO': PRO,
    'SER': SER,
    'THR': THR,
    'TRP': TRP,
    'TYR': TYR,
    'VAL': VAL,
}

PhiAtoms = {
    'ANY': ' C  '
}

PsiAtoms = {
    'ANY': ' O  '
}

OmeAtoms = {
    'ANY': ' CA '
}

ChiAtoms = {
    'ARG': (' CG ', ' CD ', ' NE ', ' CZ '),
    'ASN': (' CG ', ' OD1'),
    'ASP': (' CG ', ' OD1'),
    'CYS': (' SG ',),
    'GLN': (' CG ', ' CD ', ' OE1'),
    'GLU': (' CG ', ' CD ', ' OE1'),
    'HIS': (' CG ', ' ND1'),
    'ILE': (' CG1', ' CD1'),
    'LEU': (' CG ', ' CD1'),
    'LYS': (' CG ', ' CD ', ' CE ', ' NZ '),
    'MET': (' CG ', ' SD ', ' CE '),
    'PHE': (' CG ', ' CD1'),
    'SER': (' OG ',),
    'THR': (' OG1',),
    'TRP': (' CG ', ' CD1'),
    'TYR': (' CG ', ' CD1'),
    'VAL': (' CG1',),
}


HYDROPHOBIC = {
    ' C  ' : 1,
    ' CA ' : 1,
    ' CB ' : 1,
    ' CD ' : 1,
    ' CD1' : 1,
    ' CD2' : 1,
    ' CE ' : 1,
    ' CE1' : 1,
    ' CE2' : 1,
    ' CE3' : 1,
    ' CG ' : 1,
    ' CG1' : 1,
    ' CG2' : 1,
    ' CH2' : 1,
    ' CH3' : 1,
    ' CZ ' : 1,
    ' CZ2' : 1,
    ' CZ3' : 1,
    ' N  ' : 0,
    ' ND1' : 0,
    ' ND2' : 0,
    ' NE ' : 0,
    ' NE1' : 0,
    ' NE2' : 0,
    ' NH1' : 0,
    ' NH2' : 0,
    ' NZ ' : 0,
    ' O  ' : 0,
    ' OD1' : 0,
    ' OD2' : 0,
    ' OE1' : 0,
    ' OE2' : 0,
    ' OG ' : 0,
    ' OG1' : 0,
    ' OH ' : 0,
    ' OW1' : 0,
    ' OW2' : 0,
    ' OXT' : 0,
    ' SD ' : 0,
    ' SG ' : 0,
    ' H  ' : 1,
    ' HA ' : 1,
    ' HB ' : 1,
    ' HD1' : 1,
    ' HD2' : 1,
    ' HE ' : 1,
    ' HE1' : 1,
    ' HE2' : 1,
    ' HE3' : 1,
    ' HG ' : 1,
    ' HG1' : 1,
    ' HH ' : 1,
    ' HH2' : 1,
    ' HZ ' : 1,
    ' HZ2' : 1,
    ' HZ3' : 1,
    '1H  ' : 1,
    '1HA ' : 1,
    '1HB ' : 1,
    '1HD ' : 1,
    '1HD1' : 1,
    '1HD2' : 1,
    '1HE ' : 1,
    '1HE2' : 1,
    '1HG ' : 1,
    '1HG1' : 1,
    '1HG2' : 1,
    '1HH1' : 1,
    '1HH2' : 1,
    '1HH3' : 1,
    '1HW1' : 1,
    '1HW2' : 1,
    '1HZ ' : 1,
    '2H  ' : 1,
    '2HA ' : 1,
    '2HB ' : 1,
    '2HD ' : 1,
    '2HD1' : 1,
    '2HD2' : 1,
    '2HE ' : 1,
    '2HE2' : 1,
    '2HG ' : 1,
    '2HG1' : 1,
    '2HG2' : 1,
    '2HH1' : 1,
    '2HH2' : 1,
    '2HH3' : 1,
    '2HW1' : 1,
    '2HW2' : 1,
    '2HZ ' : 1,
    '3H  ' : 1,
    '3HB ' : 1,
    '3HD1' : 1,
    '3HD2' : 1,
    '3HE ' : 1,
    '3HG1' : 1,
    '3HG2' : 1,
    '3HH3' : 1,
    '3HZ ' : 1
    }

HYDROGEN = {
    ' H  ' : 1,
    ' HA ' : 1,
    ' HB ' : 1,
    ' HD1' : 1,
    ' HD2' : 1,
    ' HE ' : 1,
    ' HE1' : 1,
    ' HE2' : 1,
    ' HE3' : 1,
    ' HG ' : 1,
    ' HG1' : 1,
    ' HH ' : 1,
    ' HH2' : 1,
    ' HZ ' : 1,
    ' HZ2' : 1,
    ' HZ3' : 1,
    '1H  ' : 1,
    '1HA ' : 1,
    '1HB ' : 1,
    '1HD ' : 1,
    '1HD1' : 1,
    '1HD2' : 1,
    '1HE ' : 1,
    '1HE2' : 1,
    '1HG ' : 1,
    '1HG1' : 1,
    '1HG2' : 1,
    '1HH1' : 1,
    '1HH2' : 1,
    '1HH3' : 1,
    '1HW1' : 1,
    '1HW2' : 1,
    '1HZ ' : 1,
    '2H  ' : 1,
    '2HA ' : 1,
    '2HB ' : 1,
    '2HD ' : 1,
    '2HD1' : 1,
    '2HD2' : 1,
    '2HE ' : 1,
    '2HE2' : 1,
    '2HG ' : 1,
    '2HG1' : 1,
    '2HG2' : 1,
    '2HH1' : 1,
    '2HH2' : 1,
    '2HH3' : 1,
    '2HW1' : 1,
    '2HW2' : 1,
    '2HZ ' : 1,
    '3H  ' : 1,
    '3HB ' : 1,
    '3HD1' : 1,
    '3HD2' : 1,
    '3HE ' : 1,
    '3HG1' : 1,
    '3HG2' : 1,
    '3HH3' : 1,
    '3HZ ' : 1,
    ' D1 ' : 1,
    ' D2 ' : 1,
    ' D3 ' : 1
    }


def get_residue_atoms(res_name):
    res = RESIDUES[res_name]
    atoms = []
    for arec in res:
        a = Atom3d()
        a.name = arec[0]
        a.radius = arec[4]
        a.bsep0 = arec[-2]
        a.bsep1 = arec[-1]
        atoms.append(a)
    return atoms

def get_atom_with_name(atoms, start, end, name):
    for i in range(start, end):
        if atoms[i].name == name: return atoms[i]
    return None

def protein_from_pdb(p):
    """Fill in the attributes of a LinusMol

    Arguments

        o *p* - instance of LinusMol

    Result

        None
    """
    if p.filename[-3:] == '.gz':
        data = gzip.GzipFile(p.filename, 'rb').readlines()
    else:
        data = open(p.filename, 'rb').readlines()
    na = nr = -1
    a = p.atoms = []
    r = p.residue_names = []
    f = p.residue_first_atom_indices = []
    respdbnum = p.res_pdb_number = []

    nfstr = 'Unable to find atom %s in residue %i, deleting from molecule...'
    current_res = None
    aname_seen = []

    for line in data:
        if line[:4] == 'ATOM':
            na = na + 1
            rnam = line[17:20]
            rnum = line[22:27]
            if rnum <> current_res:
                alim = 0
                if len(f): alim = f[-1]
                for aptr in xrange(len(a)-1, alim, -1):
                    if a[aptr].name not in aname_seen:
#                           print (nfstr % (a[aptr].name, int(current_res)))
                        del a[aptr]
                aname_seen = []
                nr = nr + 1
                f.append(na)
                r.append(rnam)
                respdbnum.append(rnum)
                current_res = rnum
                a.extend(get_residue_atoms(rnam))
                start = na
                end = len(a)
            aname = line[12:16]
            if aname == ' OXT':
                na = na -1
                continue
            aname_seen.append(aname)
            atom = get_atom_with_name(a, start, end, aname)
            if atom is None:
                raise NameError, 'Unknown atom "%s" in resiude %s' % (aname,
                                                                      rnam)
            atom.x = float(line[30:38])
            atom.y = float(line[38:46])
            atom.z = float(line[46:54])
            atom.resnum = nr
        elif line[:4] == 'COMP':
            p.name = strip(line[10:70])
        elif line[:3] in ('TER', 'END'):
            break


    alim = 0
    if len(f): alim = f[-1]
    for aptr in xrange(len(a)-1, alim, -1):
        if a[aptr].name not in aname_seen:
#           print (nfstr % (a[aptr].name, int(current_res)))
            del a[aptr]


    p.num_atoms = na + 1
    nr = p.num_residues = nr + 1
    f.append(p.num_atoms)

    # create and initialize Zmatrix

    for i in range(nr):
        res = RESIDUES[r[i]]
        for atom_record in res:
            atom = get_atom_with_name(a, f[i], f[i+1],
                                      atom_record[0])
            if not atom: continue

            fp, fpname = atom_record[1]
            fp = fp + i

            if fp >= 0 and fp < nr:
                atom.first_parent = get_atom_with_name(a, f[fp],
                                                       f[fp+1], fpname)
                atom.fp_distance = atom.distance(atom.first_parent)
            fp, fpname = atom_record[2]
            fp = fp + i
            if fp >= 0 and fp < nr:
                atom.second_parent = get_atom_with_name(a, f[fp],
                                                       f[fp+1], fpname)
                atom.sp_angle = angle(atom, atom.first_parent,
                                             atom.second_parent)
            fp, fpname = atom_record[3]
            fp = fp + i
            if fp >= 0 and fp < nr:
                atom.third_parent = get_atom_with_name(a, f[fp],
                                                       f[fp+1], fpname)
                atom.tp_torsion = torsion(atom, atom.first_parent,
                                                 atom.second_parent,
                                                 atom.third_parent)
#           print i,atom.name,atom.fp_distance,atom.sp_angle,atom.tp_torsion

    commit_internal_coords(a)

    # setup pointers to the phi, psi, omega and chi torsions for each
    # residue

    p.phi_atoms = [None]*nr
    p.psi_atoms = [None]*nr
    p.omega_atoms = [None]*nr
    p.chi_atoms = [None]*nr


    for i in range(p.num_residues):
        name = r[i]
        try:
            phiatom = PhiAtoms[name]
        except KeyError:
            phiatom = PhiAtoms['ANY']

        try:
            psiatom = PsiAtoms[name]
        except KeyError:
            psiatom = PsiAtoms['ANY']

        try:
            omeatom = OmeAtoms[name]
        except KeyError:
            omeatom = OmeAtoms['ANY']

        chiatom = ChiAtoms.get(name, None)

        start = f[i]
        end = f[i+1]

        atom = get_atom_with_name(a, start, end, phiatom)
        if atom: p.phi_atoms[i] = atom

        atom = get_atom_with_name(a, start, end, psiatom)
        if atom: p.psi_atoms[i] = atom

        atom = get_atom_with_name(a, start, end, omeatom)
        if atom: p.omega_atoms[i] = atom

        if chiatom:
            chis = ()
            for atomname in chiatom:
                chis = chis + (get_atom_with_name(a, start, end, atomname),)
            p.chi_atoms[i] = chis

def angle(a1, a2, a3):
    """Angle in degrees between 3 atoms"""
    x1 = a1.x - a2.x
    y1 = a1.y - a2.y
    z1 = a1.z - a2.z

    x2 = a3.x - a2.x
    y2 = a3.y - a2.y
    z2 = a3.z - a2.z

    try:
        ang = (x1*x2 + y1*y2 + z1*z2)/sqrt((x1**2.0 + y1**2.0 + z1**2.0) *
                                           (x2**2.0 + y2**2.0 + z2**2.0))
    except ZeroDivisionError:
        return -999.999

    try:
        return acos(ang)*RADIANS_TO_DEGREES
    except ValueError:
        if ang < -1.0:
            return 180.0
        else:
            return 0.0

def torsion(a1, a2, a3, a4):
    """Torsion angle in degrees between 4 atoms"""
    a = a1.x - a2.x; b = a1.y - a2.y; c = a1.z - a2.z;
    d = a3.x - a2.x; e = a3.y - a2.y; f = a3.z - a2.z;
    g = a4.x - a3.x; h = a4.y - a3.y; i = a4.z - a3.z;

    ax = b*f - e*c; ay = c*d - f*a; az = a*e - b*d;
    bx = h*f - e*i; by = i*d - f*g; bz = g*e - h*d;

    try:
        ang = (ax*bx+ay*by+az*bz)/sqrt((ax**2.0 + ay**2.0 + az**2.0) *
                                       (bx**2.0 + by**2.0 + bz**2.0))
    except ZeroDivisionError:
        return -999.999
    else:
        try:
            ang = acos(ang) * RADIANS_TO_DEGREES
        except ValueError:
            if ang < -1.0:
                return 180.0
            else:
                return 0.0
        else:
            if ax*(e*bz-f*by) + ay*(f*bx-d*bz) + az*(d*by-e*bx) >= 0.0:
                return -ang
            else:
                return ang

def commit_internal_coords(atoms, start=None, end=None):
    if start is None: start = 0
    if end is None: end = len(atoms)
    for i in range(start, end):
        atoms[i].commit_internal_coords()

default_hbdpar = {
    'use_hbond': 0,
    'use_sidechain_hbond': 0,
    'hbond_distance': 3.0,
    'hbond_probe': 1.5,
    'hbond_score_short': 0.5,
    'hbond_score_long': 1.0,
    'hbond_torsion': 130.0,
    'sidechain_hbond_distance': 3.0,
    'sidechain_hbond_score': 1.0,
    'sidechain_hbond_torsion': 130.0,
    'hbond_winmin': 2,
    'hbond_winmax': 6
}

converters = {
    type(0.0): float,
    type(0): int,
    type(''): str
}

def print_hbs_chasa(fil, mol, flags, probe=1.4, ndiv=3,
                    ext_atoms=None, solv_list=None,ext_radius=1.25, title=''):
    """Prints PDB file containing:
           o number of solvating waters on backbone groups - occupancy column
           o CASA - Bfactor column
           o number of unsatisfied backbone Hbond D and A - in TER record

    Arguments

      o *filename* - filename or open file - location to which ASA
        information will be written

      o *mol* - instance of linusMol - molecule whose accessible surface
        area is to be calculated

      o *flags* - object - the object returned by the *make_asa_list*
        function in the **construct** module.

      o *probe* - float - the water probe size to use for accessible surface

      o *ndiv* - integer - a number from 1 to 5 representing the number of
        sampling points used in the area calculation.  From least (1) to
        most (5) accurate, this is either 60, 240, 960, 3840, or 15,350
        points.

      o *ext_atoms* - list of coordinates - a list of tuples, each
        tuple containing (x, y, z) of an additional atom to include
        in the ASA calculation

      o *ext_radius* - float - the radius of the atoms in ext_atoms.

      o *title* - string - An optional title for the PDB output
    """

    atoms = mol.atoms
    fai = mol.residue_first_atom_indices
    residue_names = mol.residue_names
    minres, maxres = get_res_extents(mol)


    use_ext = 0
    ext_coords = None

    use_data = 1
    data = numpy.zeros(len(atoms), 'd')

    if ext_atoms:
        ext_coords = []
        map(lambda x: map(lambda y: ext_coords.append(y), x), ext_atoms)
        ext_coords = numpy.array(ext_coords, 'd')
        use_ext = len(ext_atoms)

    tot_asa =  asa_evaluate(atoms, data, ext_coords, flags, probe,
                            ext_radius, use_data, use_ext, ndiv)

    #ff = file
    #if type(file) == type('hi'):
    ff = open(fil, 'w')

    write = ff.write
    pdbfmt = 'ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'

    p_solv_nrg = 0.0
    ap_solv_nrg = 0.0
    Gamma_p = 3.0/5.0
    Gamma_hb_oxy = 2.0
    Gamma_ap = 0.03
    CHASA = 0.0
    for i in xrange(minres,maxres):
        rname = residue_names[i]
        start = fai[i]
        end = fai[i+1]
        occ = 0.0
        for j in range(start, end):
            atom = atoms[j]
            residue_num = int(mol.res_pdb_number[atom.resnum])
            if atom.name == ' N  ':
                if solv_list[i][0][2] > 0:
                    p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i][0][2]))
                elif solv_list[i][0][2] < 0:
                    p_solv_nrg = p_solv_nrg + 1.0

                write(pdbfmt % (j+1, atom.name, rname, residue_num, atom.x, atom.y,
                            atom.z, solv_list[i][0][2],data[j]))
            elif atom.name == ' O  ':
                if solv_list[i][1][2] > 0:
                    if solv_list[i][1][3] == 0:
                        p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i][1][2]))
                    elif solv_list[i][1][3] > 0:
                        p_solv_nrg = p_solv_nrg - (Gamma_hb_oxy)
                elif solv_list[i][1][2] < 0:
                    p_solv_nrg = p_solv_nrg + 1.0

                write(pdbfmt % (j+1, atom.name, rname, residue_num, atom.x, atom.y,
                            atom.z, solv_list[i][1][2],data[j]))
            else:
                write(pdbfmt % (j+1, atom.name, rname, residue_num, atom.x, atom.y,
                            atom.z, occ,data[j]))
            if 'C' in atom.name:
                ap_solv_nrg = ap_solv_nrg + (Gamma_ap * data[j])
                CHASA = CHASA + data[j]

    aidx = j+2
    ridx = residue_num+1

    if ext_atoms:
        for x,y,z in ext_atoms:
            write(pdbfmt % (aidx, ' O  ', 'HOH', ridx, x, y, z, 0.0, 0.0))
            aidx = aidx + 1
            ridx = ridx + 1

    tot_solv_nrg = ap_solv_nrg + p_solv_nrg
    write('TER %8.3f %8.3f\n' % (CHASA,tot_solv_nrg))

    if ff <> file:
        write('END\n')
        ff.close()

    return tot_asa,tot_solv_nrg

def find_neighbors(atoms, coords, flags, probe, ext_rad,
         k, numatm, numext,nghlst):
    """
    Given an atom index k, this function calculates which atoms are close
    enough to be relevant in solvent accessible surface area calculations.
    This set of atom indices are stored in nhhlst, and the size of nghlst
    is returned.  Positive indices in nghlst represent atoms in the protein
    atom list atoms.  Negative indices (offset by one) represent atom
    coordinates stored in the external coordinate list.
    """
    numngh = 0

    a = atoms[k]
    arad = a.radius + probe + probe

    for i in range(0,k):
        if not (flags[i]): continue
        b = atoms[i]
        d = a.distance(b)
        if d < (arad + b.radius):
#           print 'found neighbor ',atoms[i].name
            nghlst.append(i)

    for i in range(k+1,numatm):
        if not (flags[i]): continue
        b = atoms[i]
        d = a.distance(b)
        if d < (arad + b.radius):
#           print 'found neighbor ',atoms[i].name
            nghlst.append(i)

    numngh = len(nghlst)
    if not numext: return numngh

    ax = a.x
    ay = a.y
    az = a.z
    for i in range(numext):
        pos = i*3
        bx = coords[pos] - ax
        by = coords[pos+1] - ay
        bz = coords[pos+2] - az
        d = sqrt(bx*bx + by*by + bz*bz)
        if d < (arad + ext_rad):
            nghlst.append(-(i+1))

    numngh = len(nghlst)
    return numngh

def tri_buried(atoms, coords, tx, ty, tz, probe, ext_rad,
                         nghlst, numngh, nghstrt):
    """
     This function determines whether one of the radial points of an atom
    is bured by steric contacts with other solvated radii.  tx, ty, tz
    are points on the solvated radii of the atom in question, i.e.
    r = atom radius + solvent radius.  This function checks whether
    any of the other solvated atoms collide with this particular point.
    If it is buried within the molecule, this function returns true.

    """
    for k in range(nghstrt,numngh):
        i = nghlst[k]
        if i >= 0:
            a = atoms[i]
            dx = a.x - tx
            dy = a.y - ty
            dz = a.z - tz
            d = dx*dx + dy*dy + dz*dz
            r = probe + a.radius
            if d < (r*r):
                nghstrt = k
                return 1,nghstrt

        else:
            pos = 3*(-(i+1))
            dx = coords[pos] - tx
            dy = coords[pos+1] - ty
            dz = coords[pos+2] - tz
            d = dx*dx + dy*dy + dz*dz
            r = ext_rad + probe
            if d < (r*r):
                nghstrt = k
                return 1,nghstrt

    for k in range(nghstrt):
        i = nghlst[k]
        if i >= 0:
            a = atoms[i]
            dx = a.x - tx
            dy = a.y - ty
            dz = a.z - tz
            d = dx*dx + dy*dy + dz*dz
            r = probe + a.radius
            if d < (r*r):
                nghstrt = k
                return 1,nghstrt

        else:
            pos = 3*(-(i+1))
            dx = coords[pos] - tx
            dy = coords[pos+1] - ty
            dz = coords[pos+2] - tz
            d = dx*dx + dy*dy + dz*dz
            r = ext_rad + probe
            if d < (r*r):
                nghstrt = k
                return 1,nghstrt

    return 0,nghstrt


def asa_evaluate(atoms, data, ext_coords, flags, probe,
                        ext_radius, use_data, use_ext, ndiv):
    """
    Given a list of atoms and a numpy array of flags, for each atom,
    this function determines the accessible surface area for each atom
    according to the flags given.  If use_data is true, then the ASA for
    each flagged atom will be stored in data (a numpy array if it used).
    If use_ext is nonzero, coords is assumed to hold a list of use_ext
    atom coordinates to include in the calculation. ext_rad is the radius
    of external water atoms in coords.  probe is the probe water radius,
    and triangles and ndiv both specify the precision of the calculation.
    To cite this algorithm, use

    Shrake, A., and J. A. Rupley.  "Environment and Exposure to Solvent of
      Protein Atoms.  Lysozyme and Insulin."  J. Mol. Bio. 79 (1973): 351-
      371.

    """
    ntrian = 960
    triarea = 4.0 * pi / ntrian
    numatm = len(atoms)
    if use_ext:
        numext = len(ext_coords)/3
    else:
        numext = 0
    PRESENT = 1
    CONTRIB = 2
    MAXNEIGHBORS = 200
    atot = 0.0
    for k in range(len(atoms)):
        nghlst = []
        if not (flags[k]):
            continue
        a = atoms[k]

        numngh = find_neighbors(atoms, ext_coords, flags, probe, ext_radius,
                 k, numatm, numext,nghlst)
        nghstrt = 0
        nts = 0
        r = probe + a.radius
        x = a.x
        y = a.y
        z = a.z
        for i in range(ntrian):
            tx = triangles[i][0]*r + x
            ty = triangles[i][1]*r + y
            tz = triangles[i][2]*r + z
            buried,nghstrt = tri_buried(atoms, ext_coords, tx, ty, tz, probe, ext_radius,
                         nghlst, numngh, nghstrt)
            if not buried:
                nts = nts + 1
        area = nts * triarea * r * r
        if use_data:
            data[k] = area

        if (flags[k] & CONTRIB):
            atot += area

    return atot

def make_hbond_list(p, wmin, wmax, hbparms):
    """
    Make list of all possible hydrogen bonding atom pairs.
    """

    fai = p.residue_first_atom_indices
    rn = p.residue_names
    atoms = p.atoms

    minres, maxres = get_res_extents(p)

    hblist = []
    hbdist = hbparms['hbond_distance']
    hbdmax = hbdist + hbparms['hbond_probe']
    hbdtor = hbparms['hbond_torsion']
    hbenes = hbparms['hbond_score_short']
    hbenel = hbparms['hbond_score_long']
    shbdis = hbparms['sidechain_hbond_distance']
    shbdmax = shbdis + hbparms['hbond_probe']
    shbtor = hbparms['sidechain_hbond_torsion']
    shbene = hbparms['sidechain_hbond_score']

    for i in xrange(minres, maxres):
        start = fai[i]
        end = fai[i+1]
        if rn[i] <> 'PRO':
            if not ((i == 0) and (minres == 0)):
                # for NH as the donor
                donor = get_atom_named(atoms, ' N  ', start, end)
                da1 = get_atom_named(atoms, ' CA ', start, end)
                da2 = get_atom_named(atoms, ' C  ', fai[i-1], start)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                                 fai[pos+1])
                            if j < 6:
                                acps.append((acp, acp1, hbdist, hbdmax, hbdtor, hbenes))
                            else:
                                acps.append((acp, acp1, hbdist, hbdmax, hbdtor, hbenel))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                                 fai[pos+1])
                            if j < 6:
                                acps.append((acp, acp1, hbdist, hbdmax, hbdtor, hbenes))
                            else:
                                acps.append((acp, acp1, hbdist, hbdmax, hbdtor, hbenel))

                # look for sidechain acceptors with wmin = 0, so GLU can Hbond itself,etc
                wmin = 0
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            scatoms = get_scacceptor(atoms, fai[pos], fai[pos+1],
                                                     rn[pos])
                            if scatoms:
                                for acp,acp1, in scatoms:
                                    acps.append((acp,acp1, shbdis, shbdmax, shbtor,
                                                 shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            scatoms = get_scacceptor(atoms, fai[pos], fai[pos+1],
                                                     rn[pos])
                            if scatoms:
                                for acp,acp1, in scatoms:
                                    acps.append((acp, acp1,shbdis, shbdmax, shbtor,
                                                 shbene))

                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ARG NH1 or NH2 as the donor
            if rn[i] == 'ARG':
                donor = get_atom_named(atoms, ' NH1', start, end)
                da1 = get_atom_named(atoms, ' CZ ', start, end)
                da2 = get_atom_named(atoms, ' NE ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

                donor = get_atom_named(atoms, ' NH2', start, end)
                da1 = get_atom_named(atoms, ' CZ ', start, end)
                da2 = get_atom_named(atoms, ' NE ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ASN OD1 as the donor(we don't know about the flip state)
            if rn[i] == 'ASN':
                donor = get_atom_named(atoms, ' OD1', start, end)
                da1 = get_atom_named(atoms, ' CG ', start, end)
                da2 = get_atom_named(atoms, ' CB ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ASN ND1 as the donor
            if rn[i] == 'ASN':
                donor = get_atom_named(atoms, ' ND2', start, end)
                da1 = get_atom_named(atoms, ' CG ', start, end)
                da2 = get_atom_named(atoms, ' CB ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ASP OD1 as the donor
            if rn[i] == 'ASN':
                donor = get_atom_named(atoms, ' OD1', start, end)
                da1 = get_atom_named(atoms, ' CG ', start, end)
                da2 = get_atom_named(atoms, ' CB ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ASP OD2 as the donor
            if rn[i] == 'ASP':
                donor = get_atom_named(atoms, ' OD2', start, end)
                da1 = get_atom_named(atoms, ' CG ', start, end)
                da2 = get_atom_named(atoms, ' CB ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for GLN OE1 as the donor
            if rn[i] == 'GLN':
                donor = get_atom_named(atoms, ' OE1', start, end)
                da1 = get_atom_named(atoms, ' CD ', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for GLN NE2 as the donor
            if rn[i] == 'GLN':
                donor = get_atom_named(atoms, ' NE2', start, end)
                da1 = get_atom_named(atoms, ' CD ', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for GLU OE1 as the donor
            if rn[i] == 'GLU':
                donor = get_atom_named(atoms, ' OE1', start, end)
                da1 = get_atom_named(atoms, ' CD ', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for GLU OE2 as the donor
            if rn[i] == 'GLU':
                donor = get_atom_named(atoms, ' OE2', start, end)
                da1 = get_atom_named(atoms, ' CD ', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for HIS NE2 as the donor
            if rn[i] == 'HIS':
                donor = get_atom_named(atoms, ' NE2', start, end)
                da1 = get_atom_named(atoms, ' CD2', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for HIS ND1 as the donor (see 5cyt, HIS 18 ND1 to PRO 30 O)
            if rn[i] == 'HIS':
                donor = get_atom_named(atoms, ' ND1', start, end)
                da1 = get_atom_named(atoms, ' CG ', start, end)
                da2 = get_atom_named(atoms, ' CB ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for LYS NZ as the donor
            if rn[i] == 'LYS':
                donor = get_atom_named(atoms, ' NZ ', start, end)
                da1 = get_atom_named(atoms, ' CE ', start, end)
                da2 = get_atom_named(atoms, ' CD ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for SER OG as the donor
            if rn[i] == 'SER':
                donor = get_atom_named(atoms, ' OG ', start, end)
                da1 = get_atom_named(atoms, ' CB ', start, end)
                da2 = get_atom_named(atoms, ' CA ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for THR OG1 as the donor
            if rn[i] == 'THR':
                donor = get_atom_named(atoms, ' OG1', start, end)
                da1 = get_atom_named(atoms, ' CB ', start, end)
                da2 = get_atom_named(atoms, ' CA ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for TYR OH as the donor
            if rn[i] == 'TYR':
                donor = get_atom_named(atoms, ' OH ', start, end)
                da1 = get_atom_named(atoms, ' CZ ', start, end)
                da2 = get_atom_named(atoms, ' CG ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

            # for ARG NH1 or NH2 as the donor
            if rn[i] == 'ARG':
                donor = get_atom_named(atoms, ' NH1', start, end)
                da1 = get_atom_named(atoms, ' CZ ', start, end)
                da2 = get_atom_named(atoms, ' NE ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))

                donor = get_atom_named(atoms, ' NH2', start, end)
                da1 = get_atom_named(atoms, ' CZ ', start, end)
                da2 = get_atom_named(atoms, ' NE ', start, end)
                acps = []
                for j in range(wmin, wmax):
                    pos = i-j
                    if pos >= minres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                 fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))

                    pos = i+j
                    if pos < maxres:
                        if hbparms['use_sidechain_hbond']:
                            acp = get_atom_named(atoms, ' O  ', fai[pos],
                                                fai[pos+1])
                            acp1 = get_atom_named(atoms, ' C  ', fai[pos],
                                             fai[pos+1])
                            acps.append((acp, acp1, shbdis, shbdmax, shbtor,shbene))
                hblist.append(((donor, da1, da2), tuple(acps)))



    return hblist

def is_hydrophobic(res, atm):
    """returns true if an atom in a residue is hydrophobic"""
    return HYDROPHOBIC[atm]

def is_hydrogen(res, atm):
    """returns true if an atom in a residue is a hydrogen"""
    return HYDROGEN.has_key(atm)

def make_asa_list(mol, hphob=1, hphil=0, hydrogens=0):
    """Construct parameters for accessible surface area scoring

    Arguments

        o *mol* - instance of linusMol - molecule whose accessible surface
          area is to be calculated

        o *hphob* - boolean - whether to include hydrophobic atoms in
          the evaluation of surface area (by default, all atoms are used
          in steric calculations except ACE, NME)

        o *hphil* - boolean - whether to include hydrophillic atoms in
          the evaluation of surface area (by default, all atoms are used
          in steric calculations except ACE, NME)

        o *hydrogens* boolean - whether to include hydrogens in the
          scoring function (as contributing to the return value or
          as steric hinderances to solvation).  If included, hydrogens
          are hydrophobic, according to *HYDROPHILLIC* in **linusRes**.

    Result

        Returns an object that can be passed to *asa_score* or
        *print_asa_score* functions in the **linusScore** module.
    """

    # first  bit: whether or not the atom is sterically present
    # second bit: whether or not the atom contributes to total ASA

    PRESENT = 1
    CONTRIB = 2

    atoms = mol.atoms
    fai = mol.residue_first_atom_indices
    flags = [0]*len(atoms)

    minres, maxres = get_res_extents(mol)

    for r in xrange(minres, maxres):
        rnam = mol.residue_names[r]
        first, last = fai[r], fai[r+1]
        for at in xrange(first, last):
            if not hydrogens and is_hydrogen(rnam, atoms[at].name):
                continue

            hydrophobic = is_hydrophobic(rnam, atoms[at].name)

            flags[at] = PRESENT

            if hphob and hydrophobic:
                flags[at] = flags[at] | CONTRIB

            if hphil and not hydrophobic:
                flags[at] = flags[at] | CONTRIB

    return numpy.array(flags, 'i')

def get_atom_named(atoms, atom_name, start, end):
    for i in range(start, end):
        atom = atoms[i]
        if atom.name == atom_name:
            return atom

def get_res_extents(p):
    """get_res_extents (linusMol p)

    Returns the residue index extents suitable for a range statement that
    scans the protein.  This information is returned as a tuple of minres,
    maxres.  Capping residues, such as ACE and NME, are ignored.
    """
    rnam = p.residue_names

    minres = 0
    maxres = p.num_residues

    if len(rnam) and rnam[minres] == 'ACE' : minres = 1
    if len(rnam) and rnam[-1] == 'NME' : maxres = maxres-1

    return minres, maxres

def angle(a1, a2, a3):
    """Angle in degrees between 3 atoms"""
    x1 = a1.x - a2.x
    y1 = a1.y - a2.y
    z1 = a1.z - a2.z

    x2 = a3.x - a2.x
    y2 = a3.y - a2.y
    z2 = a3.z - a2.z

    try:
        ang = (x1*x2 + y1*y2 + z1*z2)/sqrt((x1**2.0 + y1**2.0 + z1**2.0) *
                                           (x2**2.0 + y2**2.0 + z2**2.0))
    except ZeroDivisionError:
        return -999.999

    try:
        return acos(ang)*RADIANS_TO_DEGREES
    except ValueError:
        if ang < -1.0:
            return 180.0
        else:
            return 0.0

def norm(x1, y1, z1, x2, y2, z2, sqrt=sqrt):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    d = 1.0 / sqrt(dx*dx + dy*dy + dz*dz)
    return dx*d, dy*d, dz*d

def ztox(coordsfp, coordssp, coordstp, zmvalu, zmindx):
    """Convert internal coordinates to cartesian coordinates

    Arguments

        o *coords* - list of type 'f' and 3 elements

        o *zmvalu* - list of type 'f' and 3 elements, where each
        element is respecively the distance, angle and torsion angle of
        the atom to it's parent, second parent and third parent.

        o *zmindx* - list of type 'i' and 3 elements, where each
        element is respectively the integer index of the atoms' first,
        second and third parents.

    Returns

        atcoor - the converted x,y,z coords
    """

    fp, sp, tp = zmindx
    fpd, spa, tpt = zmvalu

    spa = spa * DEGREES_TO_RADIANS
    tpt = tpt * DEGREES_TO_RADIANS

    sina = sin(spa)
    cosa = -cos(spa)

    sint = sina * sin(tpt)
    cost = sina * cos(tpt)

    fpx, fpy, fpz = coordsfp
    spx, spy, spz = coordssp
    tpx, tpy, tpz = coordstp

    u2x, u2y, u2z = norm(fpx, fpy, fpz, spx, spy, spz)
    u1x, u1y, u1z = norm(spx, spy, spz, tpx, tpy, tpz)

    cosine = u1x*u2x + u1y*u2y + u1z*u2z
    if abs(cosine) < 1.0:
        sine = 1.0/sqrt(1.0 - cosine*cosine)
    else:
        sine = 1.0/sqrt(cosine*cosine - 1.0)

    u3x = sine * (u1y*u2z - u1z*u2y)
    u3y = sine * (u1z*u2x - u1x*u2z)
    u3z = sine * (u1x*u2y - u1y*u2x)

    u4x = cost * (u3y*u2z - u3z*u2y)
    u4y = cost * (u3z*u2x - u3x*u2z)
    u4z = cost * (u3x*u2y - u3y*u2x)

    dummy = numpy.zeros((1, 3), 'f')
    atcoor = dummy[0]

    atcoor[0] = fpx + fpd*(u2x*cosa + u3x*sint + u4x)
    atcoor[1] = fpy + fpd*(u2y*cosa + u3y*sint + u4y)
    atcoor[2] = fpz + fpd*(u2z*cosa + u3z*sint + u4z)

    return atcoor

def get_atom_named(atoms, atom_name, start, end):
    for i in range(start, end):
        atom = atoms[i]
        if atom.name == atom_name:
            return atom

def get_scacceptor(atoms, start, end, res):
    scatoms = []
    for arec in RESIDUES[res]:
        if arec[7] and arec[0] <> ' O  ':
            acp = (get_atom_named(atoms, arec[0], start, end))
            try:
                acp1 = acp.first_parent
                scatoms.append((acp,acp1))
            except AttributeError: pass
    return scatoms

def make_hbtab_loos(p, hblist, ANGLE=angle,
             ABS=abs, TORSION=torsion):
    hbtab_loos = []
    for dlist, alist in hblist:
        d1, da1, da2 = dlist
        try: Close = d1.close
        except AttributeError: continue

        for acp, acp1, hbd, hbdmax, hbtors, hbene in alist:
            # for sidechains relax orientation of donor
            if (d1.name == ' OG ' or d1.name == ' OG1' or d1.name == ' NZ ' or
                d1.name == ' OH ' or d1.name == ' NH1' or d1.name == ' NH2' or
                d1.name == ' ND1' or d1.name == ' ND2' or d1.name == ' OD2' or
                d1.name == ' SG ' or d1.name == ' NE2' or d1.name == ' OE2' or
                d1.name == ' NE1' or d1.name == ' OE1' or d1.name == ' OD1'):
                if Close(acp, hbdmax):
                    if ANGLE(d1, acp, acp1) > 90.0:
                        dist = d1.distance(acp)
#                       print dist, d1.name, p.res_pdb_number[d1.resnum], acp.name, \
#                             p.res_pdb_number[acp.resnum]
                        hbtab_loos.append((d1, acp, dist))
            elif (acp.name == ' OG ' or acp.name == ' OG1' or
                acp.name == ' OH ' or
                acp.name == ' ND1' or acp.name == ' ND2' or acp.name == ' OD2' or
                acp.name == ' SG ' or acp.name == ' NE2' or acp.name == ' OE2' or
                acp.name == ' NE1' or acp.name == ' OE1' or acp.name == ' OD1'):
                if Close(acp, hbdmax):
                    if ANGLE(d1, acp, acp1) > 90.0:
                        dist = d1.distance(acp)
#                       print dist, d1.name, p.res_pdb_number[d1.resnum], acp.name, \
#                             p.res_pdb_number[acp.resnum]
                        hbtab_loos.append((d1, acp, dist))
            else:
                if Close(acp, hbdmax):
                    # first three keep oxygen in approximate +/- 70 deg cone of NH
                    # (it's actually a square)
                    # Equivalent to theta angle >= 110 in Kortemme et al. JMB, 2003
                    # Last angle keeps N in +/- 90 deg cone of CO
                    # Equivalent to psi angle >= 90 in Kortemme
                    if ABS(TORSION(acp, d1, da1, da2)) > 130.0 and \
                    ANGLE(da1, d1, acp) > 69.0 and \
                    ANGLE(da2, d1, acp) > 69.0 and \
                    ANGLE(d1, acp, acp1) > 90.0:
                        dist = d1.distance(acp)
#                       print dist, d1.name, p.res_pdb_number[d1.resnum], acp.name, \
#                             p.res_pdb_number[acp.resnum]
                        hbtab_loos.append((d1, acp, dist))

    return hbtab_loos
def _ishbonded_loos(p,atom, hbtab_loos,numint_loos):
    """
    Tests for hydrogen bond between atom and other internal atom
    """
    rn = p.residue_names
    for rec in hbtab_loos:
        try:
            d1, acp, dist = rec
            if (atom.resnum==d1.resnum and \
                rn[atom.resnum]==rn[d1.resnum] and \
                atom.name==d1.name):
                numint_loos = numint_loos + 1
                return 1,numint_loos
            elif (atom.resnum==acp.resnum and \
                rn[atom.resnum]==rn[acp.resnum] and \
                atom.name==acp.name):
                numint_loos = numint_loos + 1
                return 1,numint_loos
        except AttributeError: continue
    return 0,numint_loos


def _close(coords1, coords2, hsd):
    """ Local version that deals with coordinates, not atom objects
    """
    x1,y1,z1 = coords1
    x2,y2,z2 = coords2
    dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    if (dist < hsd):
        return 1

    return 0

def mk_o2_solv(prot,atoms,start,end,acp,sp,tp,numvirt,
               tryall=0, watpdb=None, ext_rad=1.25):
    wat_coord = []
    #define a new single list zmatrix index
    ozmindx = [0,0,0]

    #use the O as the first parent
    ozmindx[0] = acp

    #use the C (same res) as the second parent
    ozmindx[1] = sp

    #use the CA (same res) as the third parent
    ozmindx[2] = tp

    #coords of the first parent, second and third
    coordsfp = acp.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    #define the distance, angle, and torsion for the 180 deg "on top"
    #sol1
    ozmvalu = [2.95, 180.0, 0.0]

    owat2_crds = []
    owat2_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    gotone = 0
    scale = 0.9
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
            # NH radius is 0.9, and H--O dist is 1.95
            # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat2_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx2,yy2,zz2 = owat2_crds
        format16 = 'ATOM    900  O   HOH   921    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format16 % (xx2,yy2,zz2, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx2,yy2,zz2))
            gotone = 1
    elif watpdb and tryall:
        xx2,yy2,zz2 = owat2_crds
        format16 = 'ATOM    900  O   HOH   921    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format16 % (xx2,yy2,zz2, 0.0, 0.0))

    #If the 180 worked bail out and go to NH
#   if isbump < 1 :
#       return numvirt

    #define the distance, angle, and torsion for the "lower side"
    #sol2
    ozmvalu = [2.95, 130.0, 90.0]

    owat1_crds = []
    owat1_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat1_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx1,yy1,zz1 = owat1_crds
        format15 = 'ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format15 % (xx1,yy1,zz1, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx1,yy1,zz1))
            gotone = 1
    elif watpdb and tryall:
        xx1,yy1,zz1 = owat1_crds
        format15 = 'ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format15 % (xx1,yy1,zz1, 0.0, 0.0))

    #If this worked bail out and go to NH
#   if isbump < 1 :
#       return numvirt

    #define the distance, angle, and torsion for the "upper side"
    #sol3
    ozmvalu = [2.95, 130.0, -90.0]

    owat3_crds = []
    owat3_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat3_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx3,yy3,zz3 = owat3_crds
        format17 = 'ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format17 % (xx3,yy3,zz3, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx3,yy3,zz3))
            gotone = 1
    elif watpdb and tryall:
        xx3,yy3,zz3 = owat3_crds
        format17 = 'ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx3,yy3,zz3, 0.0, 0.0))


    #If this worked bail out and go to NH
#   if isbump < 1 :
#       return numvirt

    #define the distance, angle, and torsion for the near lone pair
    #sol4
    ozmvalu = [2.95, 130.0, 0.0]

    owat4_crds = []
    owat4_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat4_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx4,yy4,zz4 = owat4_crds
        format17 = 'ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format17 % (xx4,yy4,zz4, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx4,yy4,zz4))
            gotone = 1
    elif watpdb and tryall:
        xx4,yy4,zz4 = owat4_crds
        format17 = 'ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx4,yy4,zz4, 0.0, 0.0))


    #If this worked bail out and go to NH
#   if isbump < 1 :
#       return numvirt

    #define the distance, angle, and torsion for the far lone pair
    #sol5
    ozmvalu = [2.95, -130.0, 0.0]

    owat5_crds = []
    owat5_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat5_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx5,yy5,zz5 = owat5_crds
        format17 = 'ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format17 % (xx5,yy5,zz5, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx5,yy5,zz5))
            gotone = 1
    elif watpdb and tryall:
        xx5,yy5,zz5 = owat5_crds
        format17 = 'ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx5,yy5,zz5, 0.0, 0.0))

    return numvirt,wat_coord

def mk_nh_solv(prot,atoms,start,end,donor,sp,tp,numvirt,
               tryall=0, watpdb=None, ext_rad=1.25):
    wat_coord = []
    #define a new single list zmatrix index
    nzmindx = [0,0,0]

    #use the N as the first parent
    nzmindx[0] = donor

    #use the CA as the second parent
    nzmindx[1] = sp

    #use the C (prev res) as the third parent
    nzmindx[2] = tp

    #define the distance, angle, and torsion for 180 deg
    nzmvalu = [2.95, 119.0, 180.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    gotone = 0
    scale = 0.9
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx,yy,zz))
            gotone = 1
    elif tryall and watpdb:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))

#       return numvirt

    #define the distance, angle, and torsion for CA side
    nzmvalu = [2.95, 109.0, 180.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx,yy,zz))
            gotone = 1
    elif tryall and watpdb:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))
#       return numvirt

    #define the distance, angle, and torsion for C side
    nzmvalu = [2.95, 129.0, 180.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx,yy,zz))
            gotone = 1
    elif watpdb and tryall:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))

#       return numvirt

    #define the distance, angle, and torsion for right side
    nzmvalu = [2.95, 118.5, 170.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx,yy,zz))
            gotone = 1
    elif watpdb and tryall:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))
#       return numvirt

    #define the distance, angle, and torsion for left side
    nzmvalu = [2.95, 118.5, -170.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb and (tryall or not gotone):
            watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))

        if gotone == 0:
            wat_coord.append((xx,yy,zz))
            gotone = 1
    elif watpdb and tryall:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))

#       return numvirt

    return numvirt,wat_coord

def mk_amid_solv(prot,atoms,start,end,donor,sp,tp,numvirt,
                 tryall=0,watpdb=None, ext_rad=1.25):
    wat_coord = []
    #define a new single list zmatrix index
    nzmindx = [0,0,0]

    #use the N as the first parent
    nzmindx[0] = donor

    #use the C (prev res) as the second parent
    nzmindx[1] = sp

    #use the O (prev res) as the third parent
    nzmindx[2] = tp

    #define the distance, angle, and torsion
    nzmvalu = [2.95, 180.0, 180.0]

    #coords of the first parent, second and third
    coordsfp = donor.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    nwat_crds = []
    nwat_crds = ztox(coordsfp, coordssp, coordstp, nzmvalu, nzmindx)

    #check for bumps of N hbonded virtual O with all other atoms
    isbump = 0
    scale = 0.9
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(nwat_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break
        # if no bump write out virtual water
    if not isbump:
        numvirt =  numvirt + 1
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   911    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb: watpdb.write(format14 % (xx,yy,zz, 0.0, 50.0))
        wat_coord.append((xx,yy,zz))
    elif tryall and watpdb:
        xx,yy,zz = nwat_crds
        format14 = 'ATOM    900  O   HOH   911    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format14 % (xx,yy,zz, 0.0, 0.0))


    return numvirt,wat_coord

def mk_oh_solv(prot,atoms,start,end,acp,sp,tp,numvirt,tryall=0,watpdb=None,ext_rad=1.25):
    wat_coord = []
    #define a new single list zmatrix index
    ozmindx = [0,0,0]

    #use the O as the first parent
    ozmindx[0] = acp

    #use the C (same res) as the second parent
    ozmindx[1] = sp

    #use the CA (same res) as the third parent
    ozmindx[2] = tp

    #coords of the first parent, second and third
    coordsfp = acp.coords()
    coordssp = sp.coords()
    coordstp = tp.coords()

    #define the distance, angle, and torsion for the "lower side"
    #sol2
    ozmvalu = [2.95, 130.0, 90.0]

    owat1_crds = []
    owat1_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    scale = 0.9
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat1_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx1,yy1,zz1 = owat1_crds
        format15 = 'ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb: watpdb.write(format15 % (xx1,yy1,zz1, 0.0, 50.0))
        wat_coord.append((xx1,yy1,zz1))
    elif tryall and watpdb:
        xx1,yy1,zz1 = owat1_crds
        format15 = 'ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format15 % (xx1,yy1,zz1, 0.0, 0.0))

    #If this worked bail out and go to NH
    if isbump < 1 and not tryall:
        return numvirt,wat_coord

    #define the distance, angle, and torsion for the "upper side"
    #sol3
    ozmvalu = [2.95, 130.0, -90.0]

    owat3_crds = []
    owat3_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat3_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx3,yy3,zz3 = owat3_crds
        format17 = 'ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb: watpdb.write(format17 % (xx3,yy3,zz3, 0.0, 50.0))
        wat_coord.append((xx3,yy3,zz3))
    elif tryall and watpdb:
        xx3,yy3,zz3 = owat3_crds
        format17 = 'ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx3,yy3,zz3, 0.0, 0.0))

    #If this worked bail out and go to NH
    if isbump < 1 and not tryall:
        return numvirt,wat_coord

    #define the distance, angle, and torsion for the near lone pair
    #sol4
    ozmvalu = [2.95, 130.0, 0.0]

    owat4_crds = []
    owat4_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat4_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx4,yy4,zz4 = owat4_crds
        format17 = 'ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb: watpdb.write(format17 % (xx4,yy4,zz4, 0.0, 50.0))
        wat_coord.append((xx4,yy4,zz4))
    elif tryall and watpdb:
        xx4,yy4,zz4 = owat4_crds
        format17 = 'ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx4,yy4,zz4, 0.0, 0.0))

    #If this worked bail out and go to NH
    if isbump < 1 and not tryall :
        return numvirt,wat_coord

    #define the distance, angle, and torsion for the far lone pair
    #sol5
    ozmvalu = [2.95, -130.0, 0.0]

    owat5_crds = []
    owat5_crds = ztox(coordsfp, coordssp, coordstp, ozmvalu, ozmindx)

    # check for bumps of virt water with other atoms
    isbump = 0
    for atom in atoms:
        if (prot.residue_names[atom.resnum] != 'NME') and (prot.residue_names[atom.resnum] != 'ACE'):
            atm_crds = atom.coords()
            atm_rad = scale * atom.radius
            if atom.name == ' H  ':
                # NH radius is 0.9, and H--O dist is 1.95
                # and H won't bump its own water of hydration
                hsd = 1.94
            else:
                hsd = (ext_rad + atm_rad)

            if _close(owat5_crds, atm_crds, hsd):
                #here only need the first instance of a bump
                isbump = 1
                break

    # if no bump write out virtual water
    if not isbump:
        numvirt = numvirt + 1
        xx5,yy5,zz5 = owat5_crds
        format17 = 'ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        if watpdb: watpdb.write(format17 % (xx5,yy5,zz5, 0.0, 50.0))
        wat_coord.append((xx5,yy5,zz5))
    elif tryall and watpdb:
        xx5,yy5,zz5 = owat5_crds
        format17 = 'ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
        watpdb.write(format17 % (xx5,yy5,zz5, 0.0, 0.0))

    return numvirt,wat_coord

def mk_virt_bb_cntmult_loosHbd(prot,hblist,tryall=0,watpdb=None,ext_rad=1.25):
    """
    Places oxygen in hydrogen bonding orientation to all donors and acceptors.
    Goes thru polypeptide and calls,
                mk_o2_solv, mk_nh_solv, mk_oh_solv, mk_amid_solv
    as appropriate.
    """

    numres = prot.num_residues
    atoms = prot.atoms
    rn = prot.residue_names
    fai = prot.residue_first_atom_indices

    minres,maxres = get_res_extents(prot)

    hbtab_loos = make_hbtab_loos(prot, hblist)

    wat_list = []
    cntmult_list = []
    for i in range(minres,maxres+1):
        cntmult_list.append(([i,1,0,0],[i,2,0,0]))
    numbb = 0
    numposs = 0
    bbtot = 0
    numint_loos=0
    numint_sc_loos=0
    numvirt_loos=0
    numvirt_sc_loos=0
    #backbone first
    for i in range(minres,maxres):
        start = fai[i]
        end = fai[i+1]
        acp = get_atom_named(atoms, ' O  ', start, end)
        sp = get_atom_named(atoms, ' C  ', start, end)
        tp = get_atom_named(atoms, ' CA ', start, end )

        bbtot = bbtot + 1
        numposs = numposs + 1

        hbond, numint_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_loos)

        old_numvirt_loos = numvirt_loos
        numvirt_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
             numvirt_loos,tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
        if len(wcoord) > 0:
            wat_list.append(wcoord[0])

        if hbond:
            numbb = numbb + 1
            if numvirt_loos > old_numvirt_loos:
                # if hbonded and solvent accessible just add one water
                cntmult_list[i][1][2] = 1
                # but flag this water as special for linusScore
                cntmult_list[i][1][3] = 1

        if not hbond:
            if numvirt_loos > old_numvirt_loos:
                numbb = numbb + 1
                num_solv = numvirt_loos - old_numvirt_loos
                cntmult_list[i][1][2] = num_solv
            else:
                cntmult_list[i][1][2] = -1

#               print prot.res_pdb_number[acp.resnum], prot.residue_names[acp.resnum] \
#                    , acp.name

    #construct a virtual water along the NH axis at 2.95 from the N
    for i in range(minres,maxres):
        if rn[i] <> 'PRO':
            if not ((i == 0) and (minres == 0)):
                start = fai[i]
                end = fai[i+1]
                prev = fai[i-1]
                donor = get_atom_named(atoms, ' N  ', start, end)
                sp = get_atom_named(atoms, ' CA ', start, end)
                tp = get_atom_named(atoms, ' C  ', prev, start)

                bbtot = bbtot + 1
                numposs = numposs + 1

                hbond, numint_loos = _ishbonded_loos(prot, donor, hbtab_loos,numint_loos)
                if hbond:
                    numbb = numbb + 1

                if not hbond:
                    old_numvirt_loos = numvirt_loos
                    numvirt_loos,wcoord = mk_nh_solv(prot,atoms,start,end,donor,sp,tp, \
                     numvirt_loos,tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                    if len(wcoord) > 0:
                        wat_list.append(wcoord[0])
                    if numvirt_loos > old_numvirt_loos:
                        numbb = numbb + 1
                        num_solv = numvirt_loos - old_numvirt_loos
                        cntmult_list[i][0][2] = num_solv

                    else:
                        cntmult_list[i][0][2] = -1
#                   print prot.res_pdb_number[donor.resnum], prot.residue_names[donor.resnum] \
#                     , donor.name

    # now for sidechains
    for i in range(minres,maxres):
        if rn[i] == 'ARG':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the N
            acp = get_atom_named(atoms, ' NH1', start, end)
            sp = get_atom_named(atoms, ' CZ ', start, end)
            tp = get_atom_named(atoms, ' NE ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                   numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            #construct virtual waters at 2.95 from the N
            acp = get_atom_named(atoms, ' NH2', start, end)
            sp = get_atom_named(atoms, ' CZ ', start, end)
            tp = get_atom_named(atoms, ' NE ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # try all 4 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'ASN':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OD1', start, end)
            sp = get_atom_named(atoms, ' CG ', start, end)
            tp = get_atom_named(atoms, ' CB ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)
            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                   numvirt_sc_loos, tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            donor = get_atom_named(atoms, ' ND2', start, end)
            sp = get_atom_named(atoms, ' CG ', start, end)
            tp = get_atom_named(atoms, ' CB ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, donor, hbtab_loos,numint_sc_loos)

            if not hbond:
                numvirt_sc_loos,wcoord = mk_amid_solv(prot,atoms,start,end,donor,sp,tp, \
                    numvirt_sc_loos,tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'ASP':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OD1', start, end)
            sp = get_atom_named(atoms, ' CG ', start, end)
            tp = get_atom_named(atoms, ' CB ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                    numvirt_sc_loos,tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OD2', start, end)
            sp = get_atom_named(atoms, ' CG ', start, end)
            tp = get_atom_named(atoms, ' CB ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'GLN':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OE1', start, end)
            sp = get_atom_named(atoms, ' CD ', start, end)
            tp = get_atom_named(atoms, ' CG ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                   numvirt_sc_loos, tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            donor = get_atom_named(atoms, ' NE2', start, end)
            sp = get_atom_named(atoms, ' CD ', start, end)
            tp = get_atom_named(atoms, ' CG ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, donor, hbtab_loos,numint_sc_loos)

            if not hbond:
                numvirt_sc_loos,wcoord = mk_amid_solv(prot,atoms,start,end,donor,sp,tp, \
                      numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'GLU':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OE1', start, end)
            sp = get_atom_named(atoms, ' CD ', start, end)
            tp = get_atom_named(atoms, ' CG ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                   numvirt_sc_loos, tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OE2', start, end)
            sp = get_atom_named(atoms, ' CD ', start, end)
            tp = get_atom_named(atoms, ' CG ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)

            # if C=O not hbonded try all 5 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_o2_solv(prot,atoms,start,end,acp,sp,tp, \
                  numvirt_sc_loos, tryall=tryall, watpdb=watpdb, ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'HIS':
            start = fai[i]
            end = fai[i+1]

            donor = get_atom_named(atoms, ' NE2', start, end)
            sp = get_atom_named(atoms, ' CD2', start, end)
            tp = get_atom_named(atoms, ' CG ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, donor, hbtab_loos,numint_sc_loos)
            if not hbond:
                numvirt_sc_loos,wcoord = mk_nh_solv(prot,atoms,start,end,donor,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

            donor = get_atom_named(atoms, ' ND1', start, end)
            sp = get_atom_named(atoms, ' CG ', start, end)
            tp = get_atom_named(atoms, ' CB ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, donor, hbtab_loos,numint_sc_loos)
            if not hbond:
                numvirt_sc_loos,wcoord = mk_nh_solv(prot,atoms,start,end,donor,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'LYS':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the N
            acp = get_atom_named(atoms, ' NZ ', start, end)
            sp = get_atom_named(atoms, ' CE ', start, end)
            tp = get_atom_named(atoms, ' CD ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)
            # try all 4 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'SER':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OG ', start, end)
            sp = get_atom_named(atoms, ' CB ', start, end)
            tp = get_atom_named(atoms, ' CA ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)
            # try all 4 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                   numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'THR':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OG1', start, end)
            sp = get_atom_named(atoms, ' CB ', start, end)
            tp = get_atom_named(atoms, ' CA ', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)
            # try all 4 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                     numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

        if rn[i] == 'TYR':
            start = fai[i]
            end = fai[i+1]

            #construct virtual waters at 2.95 from the O
            # one along the C=O axis
            acp = get_atom_named(atoms, ' OH ', start, end)
            sp = get_atom_named(atoms, ' CZ ', start, end)
            tp = get_atom_named(atoms, ' CE1', start, end )

            numposs = numposs + 1
            hbond, numint_sc_loos = _ishbonded_loos(prot, acp, hbtab_loos,numint_sc_loos)
            # try all 4 positions
            if not hbond:
                numvirt_sc_loos,wcoord = mk_oh_solv(prot,atoms,start,end,acp,sp,tp, \
                    numvirt_sc_loos, tryall=tryall, watpdb=watpdb,ext_rad=ext_rad)
                if len(wcoord) > 0:
                    wat_list.append(wcoord[0])

    total_loos = (numint_loos) + numvirt_loos
    bb_nonbd_loos = bbtot - numbb
    return numint_loos, numvirt_loos, bbtot, numbb,wat_list,cntmult_list
# Given a unit sphere, the array below is the coordinates of points
# lying on the sphere that can be used to estimate accessible surface area.
# There are 960 points.

triangles = [
        [ 0.142065,  0.276550,  0.950441], [-0.010435,  0.180640,  0.983494],
        [ 0.093860,  0.192683,  0.976762], [ 0.148649,  0.117947,  0.981831],
        [-0.164403,  0.075573,  0.983494], [-0.309310, -0.031464,  0.950441],
        [-0.213641, -0.017153,  0.976762], [-0.164023, -0.095418,  0.981831],
        [-0.061830,  0.090608,  0.993965], [-0.061717, -0.080548,  0.994838],
        [-0.008704,  0.012755,  0.999881], [ 0.097504,  0.028104,  0.994838],
        [-0.010177, -0.156740,  0.987587], [ 0.144419, -0.211636,  0.966620],
        [ 0.095021, -0.139247,  0.985688], [ 0.149663, -0.047667,  0.987587],
        [ 0.659784, -0.031464,  0.750796], [ 0.539741,  0.075573,  0.838432],
        [ 0.582308, -0.017153,  0.812787], [ 0.538734, -0.095418,  0.837055],
        [ 0.398310,  0.180640,  0.899288], [ 0.245163,  0.276550,  0.929201],
        [ 0.299846,  0.192683,  0.934326], [ 0.251522,  0.117947,  0.960638],
        [ 0.449659,  0.090608,  0.888593], [ 0.303643,  0.028104,  0.952371],
        [ 0.403196,  0.012755,  0.915025], [ 0.449900, -0.080548,  0.889439],
        [ 0.252866, -0.047667,  0.966326], [ 0.249395, -0.211636,  0.944993],
        [ 0.302308, -0.139247,  0.942985], [ 0.399691, -0.156740,  0.903150],
        [-0.348009,  0.509982,  0.786643], [-0.364005,  0.361769,  0.858268],
        [-0.306882,  0.449713,  0.838798], [-0.204165,  0.470843,  0.858268],
        [-0.368013,  0.203516,  0.907275], [-0.360619,  0.043725,  0.931688],
        [-0.316346,  0.133354,  0.939224], [-0.215792,  0.150880,  0.964712],
        [-0.317683,  0.294553,  0.901286], [-0.164646,  0.241277,  0.956387],
        [-0.214764,  0.314722,  0.924568], [-0.158462,  0.403205,  0.901286],
        [-0.061824,  0.255946,  0.964712], [ 0.090756,  0.351740,  0.931688],
        [-0.008845,  0.343190,  0.939224], [-0.055341,  0.416881,  0.907275],
        [-0.211240, -0.625845,  0.750796], [-0.107557, -0.536442,  0.837055],
        [-0.196321, -0.548484,  0.812787], [-0.267147, -0.475041,  0.838432],
        [ 0.000236, -0.429325,  0.903150], [ 0.106147, -0.309387,  0.944993],
        [ 0.019450, -0.332267,  0.942985], [-0.047802, -0.252840,  0.966326],
        [-0.089017, -0.448301,  0.889439], [-0.136870, -0.272499,  0.952371],
        [-0.158871, -0.370796,  0.915025], [-0.248305, -0.385678,  0.888593],
        [-0.201528, -0.191210,  0.960638], [-0.346897, -0.127467,  0.929201],
        [-0.288738, -0.208962,  0.934326], [-0.313420, -0.305039,  0.899288],
        [ 0.699503,  0.043725,  0.713291], [ 0.696646,  0.203516,  0.687943],
        [ 0.661814,  0.133354,  0.737712], [ 0.579521,  0.150880,  0.800868],
        [ 0.673594,  0.361769,  0.644511], [ 0.630591,  0.509983,  0.585041],
        [ 0.613427,  0.449713,  0.649204], [ 0.526769,  0.470843,  0.707687],
        [ 0.648047,  0.294553,  0.702334], [ 0.501791,  0.403205,  0.765266],
        [ 0.562711,  0.314722,  0.764400], [ 0.529250,  0.241277,  0.813437],
        [ 0.409434,  0.416881,  0.811526], [ 0.284882,  0.351740,  0.891696],
        [ 0.379351,  0.343190,  0.859252], [ 0.438090,  0.255946,  0.861724],
        [ 0.276003, -0.309387,  0.910001], [ 0.356751, -0.429325,  0.829704],
        [ 0.354847, -0.332267,  0.873889], [ 0.425849, -0.252840,  0.868749],
        [ 0.429643, -0.536442,  0.726386], [ 0.490790, -0.625845,  0.606170],
        [ 0.501588, -0.548484,  0.669010], [ 0.576783, -0.475041,  0.664573],
        [ 0.433318, -0.448301,  0.781832], [ 0.579301, -0.385678,  0.718097],
        [ 0.507597, -0.370796,  0.777725], [ 0.502148, -0.272499,  0.820727],
        [ 0.643341, -0.305039,  0.702184], [ 0.685915, -0.127466,  0.716431],
        [ 0.634517, -0.208961,  0.744126], [ 0.564808, -0.191210,  0.802764],
        [-0.445452,  0.498942,  0.743391], [-0.584618,  0.427737,  0.689393],
        [-0.499292,  0.427913,  0.753391], [-0.459801,  0.350916,  0.815746],
        [-0.707664,  0.342972,  0.617723], [-0.808789,  0.249821,  0.532399],
        [-0.747671,  0.259481,  0.611275], [-0.715434,  0.181892,  0.674589],
        [-0.633362,  0.349398,  0.690488], [-0.639427,  0.187484,  0.745642],
        [-0.597104,  0.271403,  0.754856], [-0.509028,  0.272874,  0.816352],
        [-0.595202,  0.107892,  0.796300], [-0.456318,  0.032883,  0.889209],
        [-0.507548,  0.111690,  0.854353], [-0.463503,  0.192697,  0.864889],
        [ 0.135742,  0.436900,  0.889209], [ 0.116528,  0.593572,  0.796300],
        [ 0.081035,  0.513335,  0.854354], [-0.010454,  0.501854,  0.864889],
        [ 0.091454,  0.732507,  0.674589], [ 0.062235,  0.844202,  0.532399],
        [ 0.030958,  0.790812,  0.611275], [-0.061373,  0.783997,  0.617723],
        [ 0.058537,  0.663770,  0.745643], [-0.094445,  0.717151,  0.690488],
        [-0.035036,  0.654954,  0.754856], [-0.068516,  0.573476,  0.816352],
        [-0.185164,  0.700322,  0.689393], [-0.302204,  0.596693,  0.743391],
        [-0.216434,  0.620933,  0.753391], [-0.159133,  0.556089,  0.815746],
        [-0.131372, -0.685178,  0.716431], [ 0.049500, -0.710273,  0.702184],
        [-0.036748, -0.667028,  0.744125], [-0.027863, -0.595645,  0.802764],
        [ 0.232067, -0.710273,  0.664573], [ 0.403843, -0.685178,  0.606170],
        [ 0.327870, -0.667028,  0.669010], [ 0.342885, -0.595645,  0.726386],
        [ 0.147936, -0.680038,  0.718097], [ 0.259470, -0.566933,  0.781832],
        [ 0.160220, -0.607844,  0.777725], [ 0.070675, -0.566933,  0.820726],
        [ 0.269715, -0.488718,  0.829704], [ 0.187471, -0.369801,  0.910001],
        [ 0.180031, -0.451560,  0.873889], [ 0.080186, -0.488718,  0.868749],
        [-0.431392, -0.137040,  0.891696], [-0.537458, -0.229271,  0.811526],
        [-0.457871, -0.228124,  0.859251], [-0.398047, -0.314627,  0.861724],
        [-0.630483, -0.318858,  0.707687], [-0.704780, -0.401265,  0.585040],
        [-0.642400, -0.407254,  0.649204], [-0.582445, -0.495343,  0.644511],
        [-0.558394, -0.320257,  0.765266], [-0.510541, -0.496059,  0.702334],
        [-0.498210, -0.409243,  0.764399], [-0.417621, -0.404861,  0.813437],
        [-0.443487, -0.574503,  0.687943], [-0.295735, -0.635419,  0.713291],
        [-0.365455, -0.567647,  0.737712], [-0.351773, -0.484629,  0.800869],
        [ 0.953363,  0.249821,  0.169376], [ 0.894196,  0.342972,  0.287722],
        [ 0.928397,  0.259481,  0.265986], [ 0.923809,  0.181892,  0.336886],
        [ 0.809496,  0.427738,  0.402190], [ 0.703004,  0.498942,  0.506796],
        [ 0.756413,  0.427913,  0.494702], [ 0.744783,  0.350916,  0.567588],
        [ 0.854704,  0.349398,  0.383930], [ 0.790242,  0.272874,  0.548687],
        [ 0.846840,  0.271403,  0.457387], [ 0.882075,  0.187484,  0.432196],
        [ 0.767608,  0.192697,  0.611266], [ 0.770620,  0.032883,  0.636446],
        [ 0.803902,  0.111690,  0.584180], [ 0.861473,  0.107892,  0.496208],
        [ 0.153263,  0.844202,  0.513647], [ 0.182623,  0.732507,  0.655807],
        [ 0.213168,  0.790812,  0.573738], [ 0.300529,  0.783997,  0.543167],
        [ 0.207697,  0.593572,  0.777518], [ 0.226770,  0.436900,  0.870456],
        [ 0.263245,  0.513335,  0.816816], [ 0.351449,  0.501854,  0.790333],
        [ 0.240944,  0.663770,  0.708065], [ 0.385598,  0.573477,  0.722799],
        [ 0.330539,  0.654954,  0.679543], [ 0.359669,  0.717152,  0.596935],
        [ 0.468598,  0.556089,  0.686426], [ 0.571420,  0.596694,  0.563414],
        [ 0.496587,  0.620933,  0.606501], [ 0.442567,  0.700323,  0.560073],
        [ 0.878628, -0.401265,  0.258841], [ 0.858857, -0.318858,  0.400866],
        [ 0.846689, -0.407254,  0.342435], [ 0.789760, -0.495343,  0.361821],
        [ 0.814449, -0.229271,  0.533018], [ 0.748706, -0.137040,  0.648582],
        [ 0.760206, -0.228124,  0.608314], [ 0.706230, -0.314627,  0.634231],
        [ 0.815396, -0.320257,  0.482250], [ 0.705126, -0.404861,  0.582139],
        [ 0.759771, -0.409243,  0.505241], [ 0.746566, -0.496059,  0.443356],
        [ 0.639671, -0.484629,  0.596620], [ 0.553582, -0.635419,  0.538322],
        [ 0.627277, -0.567647,  0.533199], [ 0.679284, -0.574503,  0.456640],
        [-0.589414,  0.770371,  0.243147], [-0.691910,  0.643752,  0.326871],
        [-0.616865,  0.717414,  0.323721], [-0.560020,  0.727500,  0.396384],
        [-0.771571,  0.492950,  0.402092], [-0.822949,  0.328280,  0.463667],
        [-0.775962,  0.416237,  0.473952], [-0.721793,  0.421259,  0.549141],
        [-0.709675,  0.576731,  0.404650], [-0.661673,  0.506272,  0.553062],
        [-0.653675,  0.584862,  0.480256], [-0.579294,  0.662218,  0.475274],
        [-0.598792,  0.506276,  0.620591], [-0.459869,  0.578830,  0.673406],
        [-0.527761,  0.585660,  0.615200], [-0.516093,  0.662828,  0.542501],
        [-0.763951, -0.314554,  0.563414], [-0.688655, -0.233611,  0.686426],
        [-0.759240, -0.236034,  0.606500], [-0.813471, -0.156790,  0.560073],
        [-0.595443, -0.144298,  0.790333], [-0.489504, -0.051880,  0.870456],
        [-0.573977, -0.057979,  0.816816], [-0.628440,  0.022998,  0.777518],
        [-0.674586, -0.149986,  0.722799], [-0.705928,  0.017632,  0.708064],
        [-0.730382, -0.069011,  0.679543], [-0.798919, -0.073461,  0.596935],
        [-0.748672,  0.096998,  0.655807], [-0.841976,  0.165058,  0.513646],
        [-0.814100,  0.089811,  0.573738], [-0.839603,  0.005978,  0.543167],
        [-0.502466,  0.829703,  0.243147], [-0.473263,  0.786703,  0.396384],
        [-0.443147,  0.835958,  0.323721], [-0.347195,  0.878983,  0.326871],
        [-0.429056,  0.722221,  0.542501], [-0.371337,  0.639244,  0.673406],
        [-0.352945,  0.704953,  0.615200], [-0.253129,  0.742154,  0.620591],
        [-0.405447,  0.780851,  0.475274], [-0.230200,  0.800706,  0.553062],
        [-0.306298,  0.821910,  0.480256], [-0.278310,  0.871092,  0.404650],
        [-0.129121,  0.825694,  0.549142], [-0.005662,  0.885991,  0.463667],
        [-0.104696,  0.874305,  0.473952], [-0.177730,  0.898184,  0.402092],
        [ 0.053324, -0.964447,  0.258841], [ 0.173326, -0.915993,  0.361821],
        [ 0.070546, -0.936889,  0.342435], [-0.016203, -0.915993,  0.400866],
        [ 0.287315, -0.841980,  0.456640], [ 0.389865, -0.747138,  0.538322],
        [ 0.299891, -0.791053,  0.533199], [ 0.218068, -0.772328,  0.596620],
        [ 0.189740, -0.876033,  0.443356], [ 0.119927, -0.804197,  0.582139],
        [ 0.104085, -0.856678,  0.505241], [ 0.000945, -0.876033,  0.482250],
        [ 0.035501, -0.772328,  0.634231], [-0.145350, -0.747137,  0.648582],
        [-0.064727, -0.791053,  0.608314], [-0.083434, -0.841980,  0.533018],
        [-0.580197, -0.796670,  0.169376], [-0.506169, -0.793915,  0.336886],
        [-0.580090, -0.769900,  0.265986], [-0.645367, -0.707614,  0.287722],
        [-0.414536, -0.762848,  0.496208], [-0.311566, -0.705594,  0.636446],
        [-0.397084, -0.707854,  0.584180], [-0.459283, -0.644525,  0.611266],
        [-0.496161, -0.753015,  0.432196], [-0.542194, -0.636371,  0.548687],
        [-0.561458, -0.689610,  0.457387], [-0.636953, -0.668498,  0.383929],
        [-0.598292, -0.565589,  0.567588], [-0.720899, -0.472719,  0.506796],
        [-0.674230, -0.548348,  0.494701], [-0.693419, -0.597841,  0.402189],
        [ 0.939203,  0.328280,  0.100644], [ 0.867672,  0.492950,  0.064390],
        [ 0.900107,  0.416238,  0.128663], [ 0.880067,  0.421259,  0.219140],
        [ 0.764766,  0.643752,  0.026779], [ 0.637523,  0.770371, -0.009616],
        [ 0.694586,  0.717414,  0.053548], [ 0.671090,  0.727501,  0.142762],
        [ 0.811827,  0.576732,  0.091203], [ 0.719976,  0.662218,  0.207610],
        [ 0.790269,  0.584862,  0.182787], [ 0.826393,  0.506272,  0.246504],
        [ 0.688492,  0.662828,  0.294343], [ 0.688586,  0.578830,  0.436812],
        [ 0.727944,  0.585660,  0.356510], [ 0.795322,  0.506276,  0.333388],
        [ 0.793692, -0.051880,  0.606104], [ 0.859336, -0.144298,  0.490632],
        [ 0.850086, -0.057979,  0.523443], [ 0.884581,  0.022998,  0.465819],
        [ 0.903889, -0.233611,  0.358344], [ 0.924434, -0.314554,  0.215588],
        [ 0.937136, -0.236034,  0.257028], [ 0.968602, -0.156790,  0.192946],
        [ 0.905342, -0.149985,  0.397316], [ 0.969804, -0.073461,  0.232558],
        [ 0.939499, -0.069011,  0.335529], [ 0.928308,  0.017632,  0.371394],
        [ 0.985924,  0.005978,  0.167088], [ 0.976435,  0.165058,  0.139033],
        [ 0.974580,  0.089811,  0.205249], [ 0.946917,  0.096998,  0.306496],
        [ 0.607263,  0.639244,  0.471804], [ 0.608542,  0.722222,  0.328744],
        [ 0.567363,  0.704953,  0.425606], [ 0.477805,  0.742154,  0.470010],
        [ 0.591397,  0.786703,  0.177053], [ 0.557655,  0.829704,  0.024750],
        [ 0.535013,  0.835958,  0.122209], [ 0.448119,  0.878983,  0.163028],
        [ 0.560284,  0.780851,  0.276323], [ 0.415585,  0.871092,  0.261700],
        [ 0.471177,  0.821910,  0.320087], [ 0.430052,  0.800706,  0.417043],
        [ 0.322184,  0.898184,  0.299105], [ 0.188464,  0.885991,  0.423675],
        [ 0.283499,  0.874305,  0.393979], [ 0.335654,  0.825695,  0.453393],
        [ 0.537751, -0.705594,  0.461477], [ 0.576908, -0.762848,  0.291960],
        [ 0.595647, -0.707854,  0.379666], [ 0.663487, -0.644525,  0.379963],
        [ 0.598108, -0.793915,  0.109393], [ 0.599900, -0.796670, -0.073737],
        [ 0.637987, -0.769900,  0.015049], [ 0.706539, -0.707614,  0.009214],
        [ 0.626586, -0.753014,  0.200898], [ 0.736836, -0.668498,  0.100914],
        [ 0.696523, -0.689610,  0.198229], [ 0.714913, -0.636371,  0.289709],
        [ 0.795922, -0.597841,  0.095369], [ 0.862509, -0.472719,  0.180596],
        [ 0.814860, -0.548347,  0.187932], [ 0.773913, -0.565589,  0.284898],
        [-0.635102,  0.754881,  0.163708], [-0.741664,  0.665917,  0.080560],
        [-0.708318,  0.686408,  0.164710], [-0.737669,  0.628238,  0.247309],
        [-0.831159,  0.556009, -0.005306], [-0.898352,  0.430169, -0.088985],
        [-0.885636,  0.464328, -0.007038], [-0.913858,  0.399457,  0.072784],
        [-0.807220,  0.584943,  0.078979], [-0.889599,  0.428997,  0.156767],
        [-0.837161,  0.522654,  0.161227], [-0.801228,  0.545692,  0.245467],
        [-0.903436,  0.359676,  0.233317], [-0.868637,  0.312790,  0.384229],
        [-0.867415,  0.385231,  0.314941], [-0.817330,  0.477436,  0.322531],
        [-0.893722,  0.147515,  0.423675], [-0.953826,  0.027443,  0.299105],
        [-0.917487,  0.054759,  0.393979], [-0.891236, -0.011527,  0.453392],
        [-0.981859, -0.096824,  0.163028], [-0.975905, -0.216788,  0.024750],
        [-0.973474, -0.193423,  0.122209], [-0.948166, -0.263883,  0.177053],
        [-0.962650, -0.069407,  0.261700], [-0.931373, -0.237045,  0.276323],
        [-0.937120, -0.139103,  0.320087], [-0.902383, -0.108539,  0.417042],
        [-0.894372, -0.303357,  0.328744], [-0.816640, -0.332417,  0.471803],
        [-0.863280, -0.271307,  0.425606], [-0.865270, -0.174351,  0.470010],
        [ 0.025419,  0.922888,  0.384229], [-0.005554,  0.972385,  0.233317],
        [-0.042483,  0.948160,  0.314941], [-0.146601,  0.935137,  0.322531],
        [-0.038797,  0.996593,  0.072784], [-0.073048,  0.993351, -0.088985],
        [-0.109492,  0.993963, -0.007038], [-0.214725,  0.976660, -0.005306],
        [-0.075148,  0.984772,  0.156767], [-0.250395,  0.964917,  0.078979],
        [-0.181476,  0.970089,  0.161227], [-0.216030,  0.945028,  0.245467],
        [-0.349695,  0.933393,  0.080560], [-0.471386,  0.866600,  0.163708],
        [-0.380933,  0.909814,  0.164710], [-0.316066,  0.915937,  0.247309],
        [ 0.125737, -0.975487,  0.180596], [ 0.266524, -0.959098,  0.095369],
        [ 0.213532, -0.958689,  0.187932], [ 0.244515, -0.926847,  0.284898],
        [ 0.401328, -0.915888,  0.009214], [ 0.523132, -0.849056, -0.073737],
        [ 0.484320, -0.874762,  0.015049], [ 0.521220, -0.846382,  0.109393],
        [ 0.353858, -0.929839,  0.100914], [ 0.472752, -0.857989,  0.200898],
        [ 0.388214, -0.899997,  0.198229], [ 0.331935, -0.897712,  0.289709],
        [ 0.500020, -0.815316,  0.291960], [ 0.460982, -0.757980,  0.461477],
        [ 0.441980, -0.812716,  0.379666], [ 0.358276, -0.852799,  0.379963],
        [-0.241049, -0.757980,  0.606104], [-0.343909, -0.815316,  0.465819],
        [-0.255929, -0.812716,  0.523443], [-0.178924, -0.852799,  0.490632],
        [-0.435542, -0.846382,  0.306496], [-0.509680, -0.849056,  0.139033],
        [-0.438935, -0.874762,  0.205249], [-0.365007, -0.915888,  0.167088],
        [-0.354854, -0.857989,  0.371394], [-0.285159, -0.929839,  0.232558],
        [-0.278254, -0.899997,  0.335529], [-0.190400, -0.897712,  0.397316],
        [-0.207127, -0.959099,  0.192946], [-0.044118, -0.975487,  0.215588],
        [-0.121865, -0.958690,  0.257028], [-0.112000, -0.926847,  0.358344],
        [-0.790032, -0.430169,  0.436811], [-0.868215, -0.399457,  0.294343],
        [-0.810741, -0.464328,  0.356510], [-0.761384, -0.556010,  0.333388],
        [-0.922091, -0.359676,  0.142762], [-0.949774, -0.312790, -0.009616],
        [-0.921265, -0.385231,  0.053548], [-0.878258, -0.477436,  0.026779],
        [-0.879125, -0.428997,  0.207610], [-0.833008, -0.545692,  0.091203],
        [-0.832720, -0.522654,  0.182787], [-0.772708, -0.584943,  0.246504],
        [-0.775352, -0.628238,  0.064390], [-0.648094, -0.754881,  0.100644],
        [-0.715744, -0.686408,  0.128662], [-0.713115, -0.665917,  0.219140],
        [ 0.790032,  0.430169, -0.436811], [ 0.761384,  0.556010, -0.333388],
        [ 0.810741,  0.464328, -0.356510], [ 0.868215,  0.399457, -0.294343],
        [ 0.713115,  0.665917, -0.219140], [ 0.648094,  0.754881, -0.100644],
        [ 0.715744,  0.686408, -0.128662], [ 0.775352,  0.628238, -0.064390],
        [ 0.772708,  0.584943, -0.246504], [ 0.833008,  0.545692, -0.091203],
        [ 0.832720,  0.522654, -0.182787], [ 0.879125,  0.428997, -0.207610],
        [ 0.878258,  0.477436, -0.026779], [ 0.949774,  0.312790,  0.009616],
        [ 0.921265,  0.385231, -0.053548], [ 0.922091,  0.359676, -0.142762],
        [ 0.906224, -0.216788, -0.362990], [ 0.966347, -0.096824, -0.238325],
        [ 0.942512, -0.193423, -0.272506], [ 0.940940, -0.263883, -0.212124],
        [ 0.994380,  0.027443, -0.102248], [ 0.988407,  0.147515,  0.035935],
        [ 0.998499,  0.054760, -0.000736], [ 0.997870, -0.011527,  0.064215],
        [ 0.987702, -0.069407, -0.140095], [ 0.993741, -0.108539,  0.026419],
        [ 0.987329, -0.139103, -0.076371], [ 0.964751, -0.237045, -0.114300],
        [ 0.980585, -0.174351,  0.089743], [ 0.936624, -0.332417,  0.110612],
        [ 0.961207, -0.271307,  0.049741], [ 0.951483, -0.303357, -0.051523],
        [ 0.031928,  0.993351, -0.110612], [ 0.064406,  0.996593,  0.051523],
        [ 0.097795,  0.993963, -0.049741], [ 0.195143,  0.976660, -0.089743],
        [ 0.097319,  0.972385,  0.212125], [ 0.128517,  0.922888,  0.362990],
        [ 0.163503,  0.948160,  0.272506], [ 0.262143,  0.935137,  0.238325],
        [ 0.130991,  0.984773,  0.114300], [ 0.295459,  0.945028,  0.140095],
        [ 0.230424,  0.970089,  0.076371], [ 0.261222,  0.964917, -0.026419],
        [ 0.388078,  0.915937,  0.102248], [ 0.497708,  0.866600, -0.035935],
        [ 0.415016,  0.909814,  0.000736], [ 0.353062,  0.933394, -0.064215],
        [ 0.868637, -0.312790, -0.384229], [ 0.903436, -0.359676, -0.233317],
        [ 0.867415, -0.385231, -0.314941], [ 0.817330, -0.477436, -0.322531],
        [ 0.913858, -0.399457, -0.072784], [ 0.898352, -0.430169,  0.088985],
        [ 0.885636, -0.464328,  0.007038], [ 0.831159, -0.556009,  0.005306],
        [ 0.889599, -0.428997, -0.156767], [ 0.807220, -0.584943, -0.078979],
        [ 0.837161, -0.522654, -0.161227], [ 0.801228, -0.545692, -0.245467],
        [ 0.741664, -0.665917, -0.080560], [ 0.635102, -0.754881, -0.163708],
        [ 0.708318, -0.686408, -0.164710], [ 0.737669, -0.628238, -0.247309],
        [-0.537751,  0.705594, -0.461477], [-0.663487,  0.644525, -0.379963],
        [-0.595647,  0.707854, -0.379666], [-0.576908,  0.762848, -0.291960],
        [-0.773913,  0.565589, -0.284898], [-0.862509,  0.472719, -0.180596],
        [-0.814860,  0.548347, -0.187932], [-0.795922,  0.597841, -0.095369],
        [-0.714913,  0.636371, -0.289709], [-0.736836,  0.668498, -0.100914],
        [-0.696523,  0.689610, -0.198229], [-0.626586,  0.753014, -0.200898],
        [-0.706539,  0.707614, -0.009214], [-0.599900,  0.796670,  0.073737],
        [-0.637987,  0.769900, -0.015049], [-0.598108,  0.793915, -0.109393],
        [-0.988407, -0.147515, -0.035935], [-0.994380, -0.027443,  0.102248],
        [-0.998499, -0.054760,  0.000736], [-0.997870,  0.011527, -0.064215],
        [-0.966347,  0.096824,  0.238325], [-0.906224,  0.216788,  0.362990],
        [-0.942512,  0.193423,  0.272506], [-0.940940,  0.263883,  0.212124],
        [-0.987702,  0.069407,  0.140095], [-0.964751,  0.237045,  0.114300],
        [-0.987329,  0.139103,  0.076371], [-0.993741,  0.108539, -0.026419],
        [-0.951483,  0.303357,  0.051523], [-0.936624,  0.332417, -0.110612],
        [-0.961207,  0.271307, -0.049741], [-0.980585,  0.174351, -0.089743],
        [-0.460982,  0.757980, -0.461477], [-0.500020,  0.815316, -0.291960],
        [-0.441980,  0.812716, -0.379666], [-0.358276,  0.852799, -0.379963],
        [-0.521220,  0.846382, -0.109393], [-0.523132,  0.849056,  0.073737],
        [-0.484320,  0.874762, -0.015049], [-0.401328,  0.915888, -0.009214],
        [-0.472752,  0.857989, -0.200898], [-0.353858,  0.929839, -0.100914],
        [-0.388214,  0.899997, -0.198229], [-0.331935,  0.897712, -0.289709],
        [-0.266524,  0.959098, -0.095369], [-0.125737,  0.975487, -0.180596],
        [-0.213532,  0.958689, -0.187932], [-0.244515,  0.926847, -0.284898],
        [-0.025419, -0.922888, -0.384229], [ 0.146601, -0.935137, -0.322531],
        [ 0.042483, -0.948160, -0.314941], [ 0.005554, -0.972385, -0.233317],
        [ 0.316066, -0.915937, -0.247309], [ 0.471386, -0.866600, -0.163708],
        [ 0.380933, -0.909814, -0.164710], [ 0.349695, -0.933393, -0.080560],
        [ 0.216030, -0.945028, -0.245467], [ 0.250395, -0.964917, -0.078979],
        [ 0.181476, -0.970089, -0.161227], [ 0.075148, -0.984772, -0.156767],
        [ 0.214725, -0.976660,  0.005306], [ 0.073048, -0.993351,  0.088985],
        [ 0.109492, -0.993963,  0.007038], [ 0.038797, -0.996593, -0.072784],
        [-0.128517, -0.922888, -0.362990], [-0.097319, -0.972385, -0.212125],
        [-0.163503, -0.948160, -0.272506], [-0.262143, -0.935137, -0.238325],
        [-0.064406, -0.996593, -0.051523], [-0.031928, -0.993351,  0.110612],
        [-0.097795, -0.993963,  0.049741], [-0.195143, -0.976660,  0.089743],
        [-0.130991, -0.984773, -0.114300], [-0.261222, -0.964917,  0.026419],
        [-0.230424, -0.970089, -0.076371], [-0.295459, -0.945028, -0.140095],
        [-0.353062, -0.933394,  0.064215], [-0.497708, -0.866600,  0.035935],
        [-0.415016, -0.909814, -0.000736], [-0.388078, -0.915937, -0.102248],
        [-0.688586, -0.578830, -0.436812], [-0.688492, -0.662828, -0.294343],
        [-0.727944, -0.585660, -0.356510], [-0.795322, -0.506276, -0.333388],
        [-0.671090, -0.727501, -0.142762], [-0.637523, -0.770371,  0.009616],
        [-0.694586, -0.717414, -0.053548], [-0.764766, -0.643752, -0.026779],
        [-0.719976, -0.662218, -0.207610], [-0.811827, -0.576732, -0.091203],
        [-0.790269, -0.584862, -0.182787], [-0.826393, -0.506272, -0.246504],
        [-0.867672, -0.492950, -0.064390], [-0.939203, -0.328280, -0.100644],
        [-0.900107, -0.416238, -0.128663], [-0.880067, -0.421259, -0.219140],
        [ 0.720899,  0.472719, -0.506796], [ 0.598292,  0.565589, -0.567588],
        [ 0.674230,  0.548348, -0.494701], [ 0.693419,  0.597841, -0.402189],
        [ 0.459283,  0.644525, -0.611266], [ 0.311566,  0.705594, -0.636446],
        [ 0.397084,  0.707854, -0.584180], [ 0.414536,  0.762848, -0.496208],
        [ 0.542194,  0.636371, -0.548687], [ 0.496161,  0.753015, -0.432196],
        [ 0.561458,  0.689610, -0.457387], [ 0.636953,  0.668498, -0.383929],
        [ 0.506169,  0.793915, -0.336886], [ 0.580197,  0.796670, -0.169376],
        [ 0.580090,  0.769900, -0.265986], [ 0.645367,  0.707614, -0.287722],
        [ 0.975905,  0.216788, -0.024750], [ 0.981859,  0.096824, -0.163028],
        [ 0.973474,  0.193423, -0.122209], [ 0.948166,  0.263883, -0.177053],
        [ 0.953826, -0.027443, -0.299105], [ 0.893722, -0.147515, -0.423675],
        [ 0.917487, -0.054759, -0.393979], [ 0.891236,  0.011527, -0.453392],
        [ 0.962650,  0.069407, -0.261700], [ 0.902383,  0.108539, -0.417042],
        [ 0.937120,  0.139103, -0.320087], [ 0.931373,  0.237045, -0.276323],
        [ 0.865270,  0.174351, -0.470010], [ 0.816640,  0.332417, -0.471803],
        [ 0.863280,  0.271307, -0.425606], [ 0.894372,  0.303357, -0.328744],
        [ 0.509680,  0.849056, -0.139033], [ 0.435542,  0.846382, -0.306496],
        [ 0.438935,  0.874762, -0.205249], [ 0.365007,  0.915888, -0.167088],
        [ 0.343909,  0.815316, -0.465819], [ 0.241049,  0.757980, -0.606104],
        [ 0.255929,  0.812716, -0.523443], [ 0.178924,  0.852799, -0.490632],
        [ 0.354854,  0.857989, -0.371394], [ 0.190400,  0.897712, -0.397316],
        [ 0.278254,  0.899997, -0.335529], [ 0.285159,  0.929839, -0.232558],
        [ 0.112000,  0.926847, -0.358344], [ 0.044118,  0.975487, -0.215588],
        [ 0.121865,  0.958690, -0.257028], [ 0.207127,  0.959099, -0.192946],
        [ 0.589414, -0.770371, -0.243147], [ 0.560020, -0.727500, -0.396384],
        [ 0.616865, -0.717414, -0.323721], [ 0.691910, -0.643752, -0.326871],
        [ 0.516093, -0.662828, -0.542501], [ 0.459869, -0.578830, -0.673406],
        [ 0.527761, -0.585660, -0.615200], [ 0.598792, -0.506276, -0.620591],
        [ 0.579294, -0.662218, -0.475274], [ 0.661673, -0.506272, -0.553062],
        [ 0.653675, -0.584862, -0.480256], [ 0.709675, -0.576731, -0.404650],
        [ 0.721793, -0.421259, -0.549141], [ 0.822949, -0.328280, -0.463667],
        [ 0.775962, -0.416237, -0.473952], [ 0.771571, -0.492950, -0.402092],
        [-0.553582,  0.635419, -0.538322], [-0.639671,  0.484629, -0.596620],
        [-0.627277,  0.567647, -0.533199], [-0.679284,  0.574503, -0.456640],
        [-0.706230,  0.314627, -0.634231], [-0.748706,  0.137040, -0.648582],
        [-0.760206,  0.228124, -0.608314], [-0.814449,  0.229271, -0.533018],
        [-0.705126,  0.404861, -0.582139], [-0.815396,  0.320257, -0.482250],
        [-0.759771,  0.409243, -0.505241], [-0.746566,  0.496059, -0.443356],
        [-0.858857,  0.318858, -0.400866], [-0.878628,  0.401265, -0.258841],
        [-0.846689,  0.407254, -0.342435], [-0.789760,  0.495343, -0.361821],
        [-0.924434,  0.314554, -0.215588], [-0.903889,  0.233611, -0.358344],
        [-0.937136,  0.236034, -0.257028], [-0.968602,  0.156790, -0.192946],
        [-0.859336,  0.144298, -0.490632], [-0.793692,  0.051880, -0.606104],
        [-0.850086,  0.057979, -0.523443], [-0.884581, -0.022998, -0.465819],
        [-0.905342,  0.149985, -0.397316], [-0.928308, -0.017632, -0.371394],
        [-0.939499,  0.069011, -0.335529], [-0.969804,  0.073461, -0.232558],
        [-0.946917, -0.096998, -0.306496], [-0.976435, -0.165058, -0.139033],
        [-0.974580, -0.089811, -0.205249], [-0.985924, -0.005978, -0.167088],
        [-0.053324,  0.964447, -0.258841], [ 0.016203,  0.915993, -0.400866],
        [-0.070546,  0.936889, -0.342435], [-0.173326,  0.915993, -0.361821],
        [ 0.083434,  0.841980, -0.533018], [ 0.145350,  0.747137, -0.648582],
        [ 0.064727,  0.791053, -0.608314], [-0.035501,  0.772328, -0.634231],
        [-0.000945,  0.876033, -0.482250], [-0.119927,  0.804197, -0.582139],
        [-0.104085,  0.856678, -0.505241], [-0.189740,  0.876033, -0.443356],
        [-0.218068,  0.772328, -0.596620], [-0.389865,  0.747138, -0.538322],
        [-0.299891,  0.791053, -0.533199], [-0.287315,  0.841980, -0.456640],
        [ 0.005662, -0.885991, -0.463667], [ 0.129121, -0.825694, -0.549142],
        [ 0.104696, -0.874305, -0.473952], [ 0.177730, -0.898184, -0.402092],
        [ 0.253129, -0.742154, -0.620591], [ 0.371337, -0.639244, -0.673406],
        [ 0.352945, -0.704953, -0.615200], [ 0.429056, -0.722221, -0.542501],
        [ 0.230200, -0.800706, -0.553062], [ 0.405447, -0.780851, -0.475274],
        [ 0.306298, -0.821910, -0.480256], [ 0.278310, -0.871092, -0.404650],
        [ 0.473263, -0.786703, -0.396384], [ 0.502466, -0.829703, -0.243147],
        [ 0.443147, -0.835958, -0.323721], [ 0.347195, -0.878983, -0.326871],
        [-0.557655, -0.829704, -0.024750], [-0.591397, -0.786703, -0.177053],
        [-0.535013, -0.835958, -0.122209], [-0.448119, -0.878983, -0.163028],
        [-0.608542, -0.722222, -0.328744], [-0.607263, -0.639244, -0.471804],
        [-0.567363, -0.704953, -0.425606], [-0.477805, -0.742154, -0.470010],
        [-0.560284, -0.780851, -0.276323], [-0.430052, -0.800706, -0.417043],
        [-0.471177, -0.821910, -0.320087], [-0.415585, -0.871092, -0.261700],
        [-0.335654, -0.825695, -0.453393], [-0.188464, -0.885991, -0.423675],
        [-0.283499, -0.874305, -0.393979], [-0.322184, -0.898184, -0.299105],
        [-0.953363, -0.249821, -0.169376], [-0.923809, -0.181892, -0.336886],
        [-0.928397, -0.259481, -0.265986], [-0.894196, -0.342972, -0.287722],
        [-0.861473, -0.107892, -0.496208], [-0.770620, -0.032883, -0.636446],
        [-0.803902, -0.111690, -0.584180], [-0.767608, -0.192697, -0.611266],
        [-0.882075, -0.187484, -0.432196], [-0.790242, -0.272874, -0.548687],
        [-0.846840, -0.271403, -0.457387], [-0.854704, -0.349398, -0.383930],
        [-0.744783, -0.350916, -0.567588], [-0.703004, -0.498942, -0.506796],
        [-0.756413, -0.427913, -0.494702], [-0.809496, -0.427738, -0.402190],
        [ 0.431392,  0.137040, -0.891696], [ 0.398047,  0.314627, -0.861724],
        [ 0.457871,  0.228124, -0.859251], [ 0.537458,  0.229271, -0.811526],
        [ 0.351773,  0.484629, -0.800869], [ 0.295735,  0.635419, -0.713291],
        [ 0.365455,  0.567647, -0.737712], [ 0.443487,  0.574503, -0.687943],
        [ 0.417621,  0.404861, -0.813437], [ 0.510541,  0.496059, -0.702334],
        [ 0.498210,  0.409243, -0.764399], [ 0.558394,  0.320257, -0.765266],
        [ 0.582445,  0.495343, -0.644511], [ 0.704780,  0.401265, -0.585040],
        [ 0.642400,  0.407254, -0.649204], [ 0.630483,  0.318858, -0.707687],
        [ 0.489504,  0.051880, -0.870456], [ 0.595443,  0.144298, -0.790333],
        [ 0.573977,  0.057979, -0.816816], [ 0.628440, -0.022998, -0.777518],
        [ 0.688655,  0.233611, -0.686426], [ 0.763951,  0.314554, -0.563414],
        [ 0.759240,  0.236034, -0.606500], [ 0.813471,  0.156790, -0.560073],
        [ 0.674586,  0.149986, -0.722799], [ 0.798919,  0.073461, -0.596935],
        [ 0.730382,  0.069011, -0.679543], [ 0.705928, -0.017632, -0.708064],
        [ 0.839603, -0.005978, -0.543167], [ 0.841976, -0.165058, -0.513646],
        [ 0.814100, -0.089811, -0.573738], [ 0.748672, -0.096998, -0.655807],
        [ 0.456318, -0.032883, -0.889209], [ 0.595202, -0.107892, -0.796300],
        [ 0.507548, -0.111690, -0.854353], [ 0.463503, -0.192697, -0.864889],
        [ 0.715434, -0.181892, -0.674589], [ 0.808789, -0.249821, -0.532399],
        [ 0.747671, -0.259481, -0.611275], [ 0.707664, -0.342972, -0.617723],
        [ 0.639427, -0.187484, -0.745642], [ 0.633362, -0.349398, -0.690488],
        [ 0.597104, -0.271403, -0.754856], [ 0.509028, -0.272874, -0.816352],
        [ 0.584618, -0.427737, -0.689393], [ 0.445452, -0.498942, -0.743391],
        [ 0.499292, -0.427913, -0.753391], [ 0.459801, -0.350916, -0.815746],
        [-0.276003,  0.309387, -0.910001], [-0.425849,  0.252840, -0.868749],
        [-0.354847,  0.332267, -0.873889], [-0.356751,  0.429325, -0.829704],
        [-0.564808,  0.191210, -0.802764], [-0.685915,  0.127466, -0.716431],
        [-0.634517,  0.208961, -0.744126], [-0.643341,  0.305039, -0.702184],
        [-0.502148,  0.272499, -0.820727], [-0.579301,  0.385678, -0.718097],
        [-0.507597,  0.370796, -0.777725], [-0.433318,  0.448301, -0.781832],
        [-0.576783,  0.475041, -0.664573], [-0.490790,  0.625845, -0.606170],
        [-0.501588,  0.548484, -0.669010], [-0.429643,  0.536442, -0.726386],
        [-0.187471,  0.369801, -0.910001], [-0.269715,  0.488718, -0.829704],
        [-0.180031,  0.451560, -0.873889], [-0.080186,  0.488718, -0.868749],
        [-0.342885,  0.595645, -0.726386], [-0.403843,  0.685178, -0.606170],
        [-0.327870,  0.667028, -0.669010], [-0.232067,  0.710273, -0.664573],
        [-0.259470,  0.566933, -0.781832], [-0.147936,  0.680038, -0.718097],
        [-0.160220,  0.607844, -0.777725], [-0.070675,  0.566933, -0.820726],
        [-0.049500,  0.710273, -0.702184], [ 0.131372,  0.685178, -0.716431],
        [ 0.036748,  0.667028, -0.744125], [ 0.027863,  0.595645, -0.802764],
        [-0.135742, -0.436900, -0.889209], [ 0.010454, -0.501854, -0.864889],
        [-0.081035, -0.513335, -0.854354], [-0.116528, -0.593572, -0.796300],
        [ 0.159133, -0.556089, -0.815746], [ 0.302204, -0.596693, -0.743391],
        [ 0.216434, -0.620933, -0.753391], [ 0.185164, -0.700322, -0.689393],
        [ 0.068516, -0.573476, -0.816352], [ 0.094445, -0.717151, -0.690488],
        [ 0.035036, -0.654954, -0.754856], [-0.058537, -0.663770, -0.745643],
        [ 0.061373, -0.783997, -0.617723], [-0.062235, -0.844202, -0.532399],
        [-0.030958, -0.790812, -0.611275], [-0.091454, -0.732507, -0.674589],
        [-0.226770, -0.436900, -0.870456], [-0.207697, -0.593572, -0.777518],
        [-0.263245, -0.513335, -0.816816], [-0.351449, -0.501854, -0.790333],
        [-0.182623, -0.732507, -0.655807], [-0.153263, -0.844202, -0.513647],
        [-0.213168, -0.790812, -0.573738], [-0.300529, -0.783997, -0.543167],
        [-0.240944, -0.663770, -0.708065], [-0.359669, -0.717152, -0.596935],
        [-0.330539, -0.654954, -0.679543], [-0.385598, -0.573477, -0.722799],
        [-0.442567, -0.700323, -0.560073], [-0.571420, -0.596694, -0.563414],
        [-0.496587, -0.620933, -0.606501], [-0.468598, -0.556089, -0.686426],
        [-0.284882, -0.351740, -0.891696], [-0.409434, -0.416881, -0.811526],
        [-0.379351, -0.343190, -0.859252], [-0.438090, -0.255946, -0.861724],
        [-0.526769, -0.470843, -0.707687], [-0.630591, -0.509983, -0.585041],
        [-0.613427, -0.449713, -0.649204], [-0.673594, -0.361769, -0.644511],
        [-0.501791, -0.403205, -0.765266], [-0.648047, -0.294553, -0.702334],
        [-0.562711, -0.314722, -0.764400], [-0.529250, -0.241277, -0.813437],
        [-0.696646, -0.203516, -0.687943], [-0.699503, -0.043725, -0.713291],
        [-0.661814, -0.133354, -0.737712], [-0.579521, -0.150880, -0.800868],
        [ 0.346897,  0.127467, -0.929201], [ 0.201528,  0.191210, -0.960638],
        [ 0.288738,  0.208962, -0.934326], [ 0.313420,  0.305039, -0.899288],
        [ 0.047802,  0.252840, -0.966326], [-0.106147,  0.309387, -0.944993],
        [-0.019450,  0.332267, -0.942985], [-0.000236,  0.429325, -0.903150],
        [ 0.136870,  0.272499, -0.952371], [ 0.089017,  0.448301, -0.889439],
        [ 0.158871,  0.370796, -0.915025], [ 0.248305,  0.385678, -0.888593],
        [ 0.107557,  0.536442, -0.837055], [ 0.211240,  0.625845, -0.750796],
        [ 0.196321,  0.548484, -0.812787], [ 0.267147,  0.475041, -0.838432],
        [ 0.348009, -0.509982, -0.786643], [ 0.204165, -0.470843, -0.858268],
        [ 0.306882, -0.449713, -0.838798], [ 0.364005, -0.361769, -0.858268],
        [ 0.055341, -0.416881, -0.907275], [-0.090756, -0.351740, -0.931688],
        [ 0.008845, -0.343190, -0.939224], [ 0.061824, -0.255946, -0.964712],
        [ 0.158462, -0.403205, -0.901286], [ 0.164646, -0.241277, -0.956387],
        [ 0.214764, -0.314722, -0.924568], [ 0.317683, -0.294553, -0.901286],
        [ 0.215792, -0.150880, -0.964712], [ 0.360619, -0.043725, -0.931688],
        [ 0.316346, -0.133354, -0.939224], [ 0.368013, -0.203516, -0.907275],
        [-0.249395,  0.211636, -0.944993], [-0.252866,  0.047667, -0.966326],
        [-0.302308,  0.139247, -0.942985], [-0.399691,  0.156740, -0.903150],
        [-0.251522, -0.117947, -0.960638], [-0.245163, -0.276550, -0.929201],
        [-0.299846, -0.192683, -0.934326], [-0.398310, -0.180640, -0.899288],
        [-0.303643, -0.028104, -0.952371], [-0.449659, -0.090608, -0.888593],
        [-0.403196, -0.012755, -0.915025], [-0.449900,  0.080548, -0.889439],
        [-0.539741, -0.075573, -0.838432], [-0.659784,  0.031464, -0.750796],
        [-0.582308,  0.017153, -0.812787], [-0.538734,  0.095418, -0.837055],
        [-0.142065, -0.276550, -0.950441], [-0.148649, -0.117947, -0.981831],
        [-0.093860, -0.192683, -0.976762], [ 0.010435, -0.180640, -0.983494],
        [-0.149663,  0.047667, -0.987587], [-0.144419,  0.211636, -0.966620],
        [-0.095021,  0.139247, -0.985688], [ 0.010177,  0.156740, -0.987587],
        [-0.097504, -0.028104, -0.994838], [ 0.061717,  0.080548, -0.994838],
        [ 0.008704, -0.012755, -0.999881], [ 0.061830, -0.090608, -0.993965],
        [ 0.164023,  0.095418, -0.981831], [ 0.309310,  0.031464, -0.950441],
        [ 0.213641,  0.017153, -0.976762], [ 0.164403, -0.075573, -0.983494]
]

def run_CHASA(file):
#   sys.stdout.write('Simulation starting time: %s\n' % time.ctime(time.time()))
    Simobj = Linus(file)
    prot = Simobj.protein
    numres = prot.num_residues
    atoms = prot.atoms
    fai = prot.residue_first_atom_indices

    hbparms = Simobj.hbdpar
    hbparms['use_hbond'] = 1
    hbparms['use_sidechain_hbond'] = 1
    wmin = 2
    wmax = numres-1
    hblist = make_hbond_list(prot, wmin, wmax, hbparms)

    numint, numvirt, bbtot, numbb,wat_list,cntmult_list = mk_virt_bb_cntmult_loosHbd(prot,hblist)
    bb_nonbd = bbtot - numbb

    sys.stdout.write('COMPND numintHbd     numSolv    num_nonHbd      total_bb_polar\n')
    sys.stdout.write('COMPND   %4d         %4d        %4d             %4d\n' % \
                    (numint, numvirt,bb_nonbd,bbtot))

    asalist  = make_asa_list(prot)
#   wat_list = None
    atot,tot_solv_nrg = print_hbs_chasa(file.split('.')[0]+'.sas', prot,
          asalist, ext_atoms=wat_list,solv_list=cntmult_list)
    int_solv_nrg = tot_solv_nrg - (2.5 * numint)
    print 'Internal Hbond + Solvation Energy = ', int_solv_nrg

#   sys.stdout.write('Simulation Finished at: %s\n' % time.ctime(time.time()))

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print USAGE
        sys.exit()

    run_CHASA(file)
