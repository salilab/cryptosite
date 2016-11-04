import numpy as np
import os, subprocess, sys
from Bio import PDB
import glob
from scipy.spatial.distance import cdist

def read_ligand_data():
    data = open('ligands.ids')
    ligands = [i.strip() for i in data.readlines()]
    data.close()

    Lxyz = {}
    Lcor = {}
    for lig in ligands:
	Lxyz[lig] = []
	data = open(lig)
	D = data.readlines()
	data.close()

	for d in D:
	    x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
	    Lxyz[lig].append(np.array([x,y,z]))
	Lxyz[lig] = np.array(Lxyz[lig])
	Lcor[lig] = D
    return Lxyz

def transform(t,xyz):

    cx, cy, cz = np.cos(t[0]),np.cos(t[1]),np.cos(t[2])
    sx, sy, sz = np.sin(t[0]),np.sin(t[1]),np.sin(t[2])

    M = np.array([np.array([cz*cy, -sy*sx*cz - sz*cx, -sy*cx*cz + sz*sx]),
                  np.array([sz*cy, -sy*sx*sz + cx*cz, -sy*cx*sz - sx*cz]),
                  np.array([sy, cy*sx, cy*cx])])

    t = np.array([t[3],t[4],t[5]])

    return np.dot(M,xyz.T).T + t



def patchmap_feature(pdb):
    '''
    Calculating PatchMap features.
    '''

    Lxyz = read_ligand_data()
    cmd = ["buildParams.pl", pdb+'.pdb', "ligands.ids", "2.0", "drug"]
    print cmd
    prc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    prc.wait()

    cmd2 = ["patch_dock.Linux", "params.txt", pdb[:-4]+".out", "7"]
    print cmd2
    prc2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE)
    prc2.wait()


    # --- Parse the data
    # - Read in protein coordinates
    data = open(pdb+'.pdb')
    D = data.readlines()
    data.close()

    XYZ = []
    Res = {}
    cnt = 0
    for d in D:
        if d[:4]=='ATOM':
            x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
            resid = int(d[22:26])
            Res[cnt] = (pdb.split('.')[0],d[17:20],resid, d[21])
            cnt += 1
            XYZ.append(np.array([x,y,z]))
    XYZ = np.array(XYZ)


    # - Create a dictionary with residue contact frequencies
    Counts = dict([ (c,0) for c in set(Res.values()) ])

    for i,lig in enumerate(ligands):

        try: data = open(pdb+'.pdb'+str(i)+'.res')
        except IOError: continue
        D = data.readlines()
        data.close()

        output = open('tr','w')
        for d in D[3:-2]:
            d = d.strip().split('|')
            output.write(d[0]+' '+d[-1]+'\n')
        output.close()


        proc = subprocess.Popen(['ligand_score_multiple',\
                                 pdb+'.pdb', 'Ligands/%s.mol2' % lig.split('/')[-1].split('.')[0], 'tr'],stdout=subprocess.PIPE)

        proc.wait()



        print i,lig,pdb
        data = open('mol2_score.res')#pdb+str(i)+'.res')
        LS = [float(l.strip().split()[-1]) for l in data.readlines()]
        data.close()

        data = open(pdb+'.pdb'+str(i)+'.res')
        D = data.readlines()
        data.close()

        h = []
        for a,d in enumerate(D[3:-2]):
            d = d.strip().split()
            t = np.array([float(j) for j in d[-6:]])
            xyz = transform(t, Lxyz[lig])

            ligsc = LS[a] #ligand_score(lig,xyz,Lcor[lig],pdb)
            if ligsc>0.: continue

            Dists = cdist(XYZ, xyz)
            Dists = set(np.argwhere(Dists<3.5)[:,0])
            Dists = [Res[j] for j in Dists]
            for j in Dists: Counts[j] += 1


    output = open(pdb+'.pdb.ptm','w')
    for c in Counts: output.write('\t'.join([str(j) for j in c])+'\t'+str(Counts[c])+'\n')
    output.close()

if __name__=='__main__':

    patchmap_feature('XXXX') #to run test.pdb
