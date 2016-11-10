#!/usr/bin/env python

"""Add features relating to pockets detected by fpocket."""

from __future__ import print_function, absolute_import
import sys, os, glob
from numpy import linalg, array, shape, argwhere
from scipy import spatial
import subprocess
import optparse
from operator import itemgetter

def get_cnc(pfil,ifil):

    data = open(pfil)
    D = data.readlines()
    data.close()

    Pockets = {}
    PocketC = []
    Protins = {}
    ProtinC = []
    pt,pc = 0,0
    for d in D:
        if d[17:20]=='STP' and d[:6]=='HETATM':
            pn = int(d[22:26])
            Pockets[pc] = pn
            PocketC.append(array([float(d[30:38]), float(d[38:46]), float(d[46:54])]))
            pc += 1
        elif d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            if cid == ' ':
                cid='A'
            coords = [float(d[30:38]), float(d[38:46]), float(d[46:54])]
            ProtinC.append(array(coords))
            Protins[pt] = (atom,res,resid,cid)
            pt+=1

    # if no pockets have been found
    if PocketC==[]:
        Residues = {}
        for r in set([i[1:] for i in Protins.values()]):
            Residues[r] = 0.
        return Residues

    PocketC = array(PocketC)
    ProtinC = array(ProtinC)


    PocketInfo = {}
    data = open(ifil)
    Z = data.read().split('\n\n')
    data.close()

    for p in Z:
        p = p.split('\n')
        pn = ''
        for i in p:
            if i[:6]=='Pocket': pn = int(i.split(':')[0][6:])
            if '\tDruggability Score :' in i: PocketInfo[pn] = float(i.split(':')[1])

    Dists = spatial.distance.cdist(PocketC, ProtinC)
    Sel = argwhere(Dists<=5.)

    Residues = {}
    for p in Sel:
        pc,pt = p[0],p[1]
        pci,pti = Pockets[pc], Protins[pt][1:]

        if pci not in PocketInfo: PocketInfo[pci] = 'NAN'
        if pti not in Residues: Residues[pti] = PocketInfo[pci]
        else: Residues[pti] = max([Residues[pti], PocketInfo[pci]])

    for r in set([i[1:] for i in Protins.values()]):
        if r not in Residues: Residues[r] = 0.

    return Residues



def run_fpocket(Snap):
    subprocess.check_call(["fpocket", "-f", Snap])

def pocket_features():
    RES = {}

    out = open('pockets.out','w')
    snaps = []

    data = open('SnapList.txt')
    Data = data.readlines()
    soap_scores = [float(i.strip().split()[-1]) for i in Data]
    DRS = [i.strip().split()[0] for i in Data if float(i.strip().split()[-1]) < min(soap_scores)*.4]
    data.close()
    x=0

    print('DRS: ', DRS)

    for dr in DRS:
        run_fpocket(dr)
        fils = glob.glob(dr.rsplit('.',1)[0]+'_out/*_out.pdb')
        print(fils)
        if len(fils)==0: continue
        fil = fils[0]
        pdbfil = fil
        inffil = fil.rsplit('_',1)[0]+'_info.txt'
        print(pdbfil, inffil)
        res = get_cnc(pdbfil,inffil)
        snaps.append(fil)

        if x==0:
            for r in res: RES[r] = [res[r]]
            x += 1
        else:
            for r in res:
                if r in RES: RES[r].append(res[r])
                else: RES[r] = [res[r]]

        x += 1

    out.write('\t'.join(['Res','ResID','ChainID']+snaps)+'\n')

    # --- r[2] chainID added to pockets.out file. Eventualy columns adressed later have to be changed.
    for r in sorted(RES.keys(), key=itemgetter(1)):
        out.write('\t'.join([str(i) for i in [r[0],r[1],r[2]]+RES[r]])+'\n')
    out.close()

def parse_args():
    usage = """%prog [opts]

Add features relating to pockets detected by fpocket. This tool should be
run in the same directory that 'cryptosite soap' was previously run in.
For each AllosMod model previously scored well by SOAP, fpocket is used
to detect pockets and features derived from this output are stored in
a pockets.out file.
"""
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

def main():
    parse_args()
    pocket_features()

if __name__ == '__main__':
    main()
