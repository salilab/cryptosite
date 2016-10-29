#!/usr/bin/env python

from modeller import *
import os, sys, warnings
from numpy import dot, transpose, linalg, sqrt, array
import numpy
from operator import itemgetter
import subprocess
from CHASA import *

warnings.filterwarnings("ignore")


def get_sas(pdb,probe):

    # Read the PDB file
    env = environ()
    mdl = model(env)
    mdl.read(file=pdb)

    # Calculate atomic accessibilities (in Biso) with appropriate probe_radius
    myedat = energy_data()
    myedat.radii_factor = 1.6
    mdl.write_data(edat=myedat, output='PSA ATOMIC_SOL',
                   psa_integration_step=0.05, probe_radius=probe)

    mdl.write(file=pdb.rsplit('.',1)[0]+'.sas')

    # read SAS
    data = open('%s.sas' % (pdb.rsplit('.',1)[0], ))
    D = data.readlines()
    data.close()

    Sas = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            if cid == ' ':
                cid='A'
            Sas[(atom,res,resid,cid)] = float(d[60:66])
    return Sas


def get_sas2(pdb):
    '''
    Calculate accesible surface area using CHASA algorithm.
    '''
    #TODO: add option for calculating that in Modeller with varying rolling sphere radius.

    # do SAS
    run_CHASA(pdb)

    # read SAS
    data = open('%s.sas' % (pdb.rsplit('.',1)[0], ))
    D = data.readlines()
    data.close()

    Sas = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            if cid == ' ':
                cid='A'
            Sas[(atom,res,resid,cid)] = float(d[60:66])
    return Sas


def protein_protrusion(pdb):

    data = open(pdb)
    D = data.readlines()
    data.close()

    XYZ = []
    for d in D:
        if 'ATOM'==d[:4]:
            x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
            XYZ.append(numpy.array([x,y,z]))
    XYZ = numpy.array(XYZ)


    out = open(pdb+'.prt','w')
    for d in D:
        if 'ATOM'==d[:4]:

            xyz = numpy.array([float(d[30:38]), float(d[38:46]), float(d[46:54])])
            distances = numpy.sqrt(numpy.sum((XYZ-xyz)**2,axis=1))
            mins = set(numpy.argwhere(8. <= distances).transpose()[0])
            maxs = set(numpy.argwhere(distances <= 12.).transpose()[0])
            natoms = str(len(mins&maxs))+'.00'
            if len(natoms)>6: natoms = '999.00'
            elif len(natoms)<6: natoms = (6-len(natoms))*' ' + natoms

            out.write(d[:60]+str(natoms)+d[66:])
        else: out.write(d)
    out.close()



def get_prt(apo):
    '''
    Calculate protein protrusion using Andreas Fisher equation.
    '''

    # do protrusion
    protein_protrusion(apo)

    # read protrusion
    data = open('%s.prt' % (apo))
    D = data.readlines()
    data.close()

    Prt = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            if cid == ' ':
                cid='A'
            Prt[(atom,res,resid,cid)] = float(d[60:66])
    return Prt


def protein_convexity(pdb):
    '''
    Calculate protein convexity as in Andras Fisher paper.
    '''

    data = open(pdb.rsplit('.',1)[0]+'.sas')
    D = data.readlines()
    data.close()

    XYZ = {}
    for d in D:
        if 'ATOM'==d[:4] and d[17:20]!='HOH':
            x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
            atid,rsid,sas,cid = int(d[6:11]), int(d[22:26]), float(d[60:66]), d[21]
            if cid == ' ':
                cid='A'
            if (rsid,cid) not in XYZ: XYZ[rsid,cid] = [[atid,x,y,z,sas]]
            else: XYZ[(rsid,cid)].append([atid,x,y,z,sas])

    SurfRes = {}
    for res in XYZ:
        sass = [i[-1] for i in XYZ[res]]
        if max(sass)>2.: SurfRes[res] = XYZ[res]

    Centroids = {}
    for res in SurfRes:
        C = SurfRes[res]
        Cs = numpy.array([i for i in C if i[-1]>=2.])
        aac = numpy.mean(C, axis=0)[1:-1]
        sac = numpy.mean(Cs, axis=0)[1:-1]
        slv = 2*sac - aac
        Centroids[res] = (sac,slv)

    SurfNet = {}
    RX = SurfRes.keys()
    for r1 in xrange(len(RX)-1):
        F1 = SurfRes[RX[r1]]
        C1 = numpy.array([numpy.array([c[1],c[2],c[3]]) for c in F1 if c[4]>2.])
        for r2 in xrange(r1+1,len(RX)):
            F2 = SurfRes[RX[r2]]
            C2 = numpy.array([numpy.array([c[1],c[2],c[3]]) for c in F2 if c[4]>2.])

            dx = numpy.subtract.outer(C1[:,0], C2[:,0])
            dy = numpy.subtract.outer(C1[:,1], C2[:,1])
            dz = numpy.subtract.outer(C1[:,2], C2[:,2])

            distance = numpy.sqrt(dx**2 + dy**2 + dz**2)
            if numpy.min(distance) < 6.:
                a1,a2 = RX[r1], RX[r2]
                dslv = numpy.linalg.norm(Centroids[a1][1]-Centroids[a2][1])
                dexp = numpy.linalg.norm(Centroids[a1][0]-Centroids[a2][0])

                fc = (dslv/dexp)-1

                if a1 not in SurfNet: SurfNet[a1] = [fc]
                else: SurfNet[a1].append(fc)
                if a2 not in SurfNet: SurfNet[a2] = [fc]
                else: SurfNet[a2].append(fc)

    out = open(pdb+'.con','w')
    for d in D:
        if 'ATOM'==d[:4]:
            res,cid = int(d[22:26]),d[21]
            if cid == ' ':
                cid='A'
            if (res,cid) in SurfNet:
                cnvx = str(round(100*numpy.mean(SurfNet[(res,cid)]),2))
                if len(cnvx)>6: cnvx = '999.00'
                elif len(cnvx)<6: cnvx = (6-len(cnvx))*' ' + cnvx
            else: cnvx = '  0.00'

            out.write(d[:60]+str(cnvx)+d[66:])
        else: out.write(d)
    out.close()



def get_cvx(apo):
    '''
    Calculate protein convexity as in Andras Fisher paper.
    '''

    # do convexities
    protein_convexity(apo)

    # read convexities
    data = open('%s.con' % (apo))
    D = data.readlines()
    data.close()

    Cvx = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            if cid == ' ':
                cid='A'
            Cvx[(atom,res,resid,cid)] = float(d[60:66])
    return Cvx



def main():
    out = open('am_features.out','w')
    #out = open('am_features_%s.out' % sys.argv[-1],'w')

    snaps = []

    data = open('SnapList.txt')
    Data = data.readlines()
    soap_scores = [float(i.strip().split()[-1]) for i in Data]
    DRS = [i.strip().split()[0] for i in Data if float(i.strip().split()[-1]) < min(soap_scores)*.4]
    data.close()
    x=0

    print 'DRS: ', DRS

    for dr in DRS:
	   # if 'pm.pdb' not in dr: continue

	'''
	Gather bioinformatics features (no neighborhood yet).
	'''

	pdb = dr #sys.argv[-1]
	chain = pdb[-1]
	pdb = pdb
	print pdb

	try:
	    sasa14 = get_sas(pdb,probe=1.4)
	    sasa30 = get_sas(pdb,probe=3.)
	    prta = get_prt(pdb)
	    cvxa = get_cvx(pdb)
	except: continue


	#outf = open(pdb+'.feat','w')

	# read pdb fragment PDB
	data = open(pdb)
	D = data.readlines()
	data.close()

	RES = {}
	for d in D:
	    d = d.strip()
	    if d[:4]=='ATOM':
		atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
		if cid == ' ':
		    cid='A'
		if atom==' OXT': continue

		p = (d[17:20], int(d[22:26]), d[21])
		if d[21]== ' ':
		    p = (d[17:20], int(d[22:26]), 'A')
		if p not in RES: RES[p] = {'sas14':[],'sas30':[],'prt':[],'cvx':[]}
		sasa14i = sasa14[(atom,res,resid,cid)]
		sasa30i = sasa30[(atom,res,resid,cid)]
		prtai = prta[(atom,res,resid,cid)]
		cvxai = cvxa[(atom,res,resid,cid)]

		RES[p]['sas14'].append(sasa14i)
		RES[p]['sas30'].append(sasa30i)
		RES[p]['prt'].append(prtai)
		RES[p]['cvx'].append(cvxai)


	for p in RES:
	    L = [str(i) for i in list(p)+[numpy.mean(RES[p]['sas14']),numpy.mean(RES[p]['sas30']),numpy.mean(RES[p]['prt']),numpy.mean(RES[p]['cvx'])]]
	    out.write( '\t'.join(L)+'\n' )
    out.close()

if __name__ == '__main__':
    main()
