from __future__ import print_function, absolute_import
import os, subprocess, sys
from Bio import PDB



'''
Wimley-White whole residue hydrophobicity interface scale
'''
hydrophobicity = {'A':0.17,
'R':0.81,
'N':0.42,
'D':1.23,
'C':-0.24,
'Q':0.58,
'E':2.02,
'G':0.01,
'H':0.96,
'I':-0.31,
'L':-0.56,
'K':0.99,
'M':-0.23,
'F':-1.13,
'P':0.45,
'S':0.13,
'T':0.14,
'W':-1.85,
'Y':-0.94,
'V':0.07,
'X':0.}


charge = {'A':0,
'R':1,
'N':0,
'D':-1,
'C':0,
'Q':0,
'E':-1,
'G':0,
'H':0.5,
'I':0,
'L':0,
'K':1,
'M':0,
'F':0,
'P':0,
'S':0,
'T':0,
'W':0,
'Y':0,
'V':0,
'X':0}



def HydChrSSE(pdb,chain):
    '''
    Predict hydrophobicity, charges, and SSE.
    '''
    command = ["mkdssp", "-i", pdb+".pdb", "-o", "%s.ssp" % pdb]
    prc = subprocess.Popen(command, stdout=subprocess.PIPE)
    prc.wait()

    # --- added for loop, change of DIC for cid 3to1 residue identification
    data = open("%s.pdb" % pdb)
    D = data.readlines()
    data.close()

    DIC = {}
    for d in D:
        if d[:4]=='ATOM' and d[21]==chain:
            res,resid = d[17:20],int(d[22:26])
            DIC[str(resid)+chain] = res


    data = open("%s.ssp" % pdb)
    D = data.read().split(' #  RESIDUE')
    data.close()

    # --- added cid -> chainID, fields change in XXX.hcs. XXX.ssp contains all chains but DIC only the ones from query - second if clausel "cid in chains".
    out = open("%s_mdl%s.hcs" % (pdb[:3],chain), "w")
    D = D[-1].split('\n')[1:-1]
    for d in D:
        rid, cid, res, ss = d[5:11],d[11],d[13],d[16]
        if res=='!': continue
        if ss==' ': ss='-'
        if cid == chain:
            try: L = '\t'.join([ rid, res, ss, str(hydrophobicity[res]), str(charge[res]) ]) + '\n'
            except KeyError:  L = '\t'.join([ rid, res, ss, \
                                      str(hydrophobicity[PDB.Polypeptide.three_to_one(DIC[str(int(rid))+cid])]),\
                                      str(charge[PDB.Polypeptide.three_to_one(DIC[str(int(rid))+cid])]) ]) + '\n'
            out.write(L)
    out.close()

    # os.system("rm "%s.ssp" % pdb)
