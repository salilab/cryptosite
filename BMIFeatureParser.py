from SiteCrypt import PATH2BIOPYTHON, PATH2FPOCKET, PATH2CONCAVITY

import os, sys, warnings
sys.path.append(PATH2BIOPYTHON)
from Bio.PDB.PDBParser import PDBParser
from Bio import PDB
from Bio import SeqIO
from numpy import dot, transpose, linalg, sqrt, array
import numpy
from operator import itemgetter
import subprocess
from CHASA import *
from modeller import *


warnings.filterwarnings("ignore")

parser = PDBParser()



def get_sas(pdb,prog='CHASA',probe=1.4):
    '''
    Calculate accesible surface area using CHASA algorithm.
    '''
    #TODO: add option for calculating that in Modeller with varying rolling sphere radius.

    # do SAS
    if prog == 'CHASA': run_CHASA(pdb+'.pdb')
    elif prog =='MOD':
        # Read the PDB file
        env = environ()
        mdl = model(env)
        mdl.read(file=pdb+'.pdb')

        # Calculate atomic accessibilities (in Biso) with appropriate probe_radius
        myedat = energy_data()
        myedat.radii_factor = 1.6
        mdl.write_data(edat=myedat, output='PSA ATOMIC_SOL',
                       psa_integration_step=0.05, probe_radius=probe)

        mdl.write(file=pdb+'.sas')

    # read SAS
    data = open('%s.sas' % (pdb))
    D = data.readlines()
    data.close()

    Sas = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
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
    protein_protrusion(apo+'.pdb')

    # read protrusion
    data = open('%s.pdb.prt' % (apo))
    D = data.readlines()
    data.close()

    Prt = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            Prt[(atom,res,resid,cid)] = float(d[60:66])
    return Prt


def get_cnc(apo):
    '''
    Find pockets using Fpocket algorithm.
    '''

    global PATH2FPOCKET

    command = [PATH2FPOCKET+"fpocket", "-f", apo+'.pdb']
    prc = subprocess.Popen(command, stdout=subprocess.PIPE)
    prc.wait()

    data = open('%s_out/%s_out.pdb' % (apo,apo))
    D = data.readlines()
    data.close()

    Pockets = {}
    for d in D:
        if d[17:20]=='STP' and d[:6]=='HETATM':
            pn = int(d[22:26])
            if pn not in Pockets: Pockets[pn] = [(float(d[30:38]), float(d[38:46]), float(d[46:54]))]
            else: Pockets[pn].append( (float(d[30:38]), float(d[38:46]), float(d[46:54])) )

    PocketInfo = {}
    data = open('%s_out/%s_info.txt' % (apo,apo))
    Z = data.read().split('\n\n')
    data.close()


    for p in Z:
        p = p.split('\n')
        pn = ''
        for i in p:
            if i[:6]=='Pocket': pn = int(i.split(':')[0][6:])
            if 'Druggability Score' in i.strip().split(':')[0].strip():
                PocketInfo[pn] = float(i.split(':')[1])


    Atoms = {}
    Residues = {}
    for d in D:
        if d[:4]=='ATOM':
            atom = d[6:16]
            coords = (float(d[30:38]), float(d[38:46]), float(d[46:54]))
            adist = {}
            for p in Pockets:
                mini=1000.
                for pi in Pockets[p]:
                    dist = linalg.norm( array(coords)-array(pi) )
                    if dist<mini: mini=dist
                if mini<=5.: adist[p] = mini
            atom, res, resid, cid  = d[12:16], d[17:20], int(d[22:26]), d[21]
            Residues[(res,resid,cid)] = 0.

            if len(adist)>0:
                if max([PocketInfo[i] for i in adist])>Residues[(res,resid,cid)]: Residues[(res,resid,cid)]=max([PocketInfo[i] for i in adist])
    print Residues
    return Residues, ('1','1')




def get_cnc2(apo):
    '''
    Find pockets using Concavity algorithm.
    '''

    global PATH2CONCAVITY

    # do concavities
    cmd = [PATH2CONCAVITY+'concavity', '-grid_method', 'pocketfinder', '-extraction_method', 'none', '-print_grid_pdb', '1', apo+'.pdb', 'cnc']
    prc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    prc.wait()

    # read concavities
    data = open('%s_cnc_residue.pdb' % (apo))
    D = data.readlines()
    data.close()

    Cnc = {}
    for d in D:
        d = d.strip()
        if d[:4]=='ATOM':
            atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
            Cnc[(atom,res,resid,cid)] = float(d[60:66])

    # extract grid points
    data = open('%s_BS.pdb' % (apo))
    D = data.readlines()
    data.close()

    XYZ = []
    for d in D:
        if 'ATOM'==d[:4]:
            x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
            XYZ.append(numpy.array([x,y,z]))
    XYZ = numpy.array(XYZ)


    data = open('%sp_cnc_pocket.pdb' % (apo))
    D = data.readlines()
    data.close()
    Points = []
    for d in D:
        if 'HETATM'==d[:6]:
            hind = d.index(' H ')
            xyz = numpy.array([float(d[19+hind:27+hind]), float(d[27+hind:35+hind]), float(d[35+hind:43+hind])])

            distances = numpy.sqrt(numpy.sum((XYZ-xyz)**2,axis=1))
            if numpy.min(distances) <= 4. and float(d[49+hind:55+hind])>=3.:
                Points.append(float(d[49+hind:55+hind]))

    print
    print 'HISTOGRAM %s:' % apo,achain
    print numpy.histogram(Points, bins=10)
    print


    return Cnc, (str(len(Points)), str(numpy.mean(Points)), str(numpy.std(Points)))



def protein_convexity(pdb):
    '''
    Calculate protein convexity as in Andras Fisher paper.
    '''

    data = open(pdb+'.sas')
    D = data.readlines()
    data.close()

    XYZ = {}
    for d in D:
        if 'ATOM'==d[:4] and d[17:20]!='HOH':
            x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
            atid,rsid,sas,cid = int(d[6:11]), int(d[22:26]), float(d[60:66]), d[21]
            if (rsid,cid) not in XYZ: XYZ[(rsid,cid)] = [[atid,x,y,z,sas]]
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
            if cid==' ':
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
            Cvx[(atom,res,resid,cid)] = float(d[60:66])
    return Cvx


def get_hcs(apo,achain):
    '''
    Retrieve hydrophobicity, charge, and SSEs.
    '''

    # ---  change of DIC for cid 3to1 residue identification
    data = open("%s.pdb" % apo)
    D = data.readlines()
    data.close()

    DIC = {}
    for d in D:
        if d[:4]=='ATOM' and d[21]==achain:
            res,resid = d[17:20],int(d[22:27])
            DIC[str(resid)+achain] = res



    # read hydrophobicity, charge, SSE
    data = open('%s%s.hcs' % (apo, achain))
    D = data.readlines()
    data.close()

    Hcs = {}
    for d in D:
        d = d.strip().split()
        if len(d)>0:
            if d[1]=='X': res, resid = 'UNK', int(d[0])
            else:
                try: res, resid = PDB.Polypeptide.one_to_three(d[1]), d[0]
                except KeyError:
                    res, resid = DIC[str(int(d[0]))+achain], d[0]
            Hcs[(res,resid)] = (d[2],float(d[3]),float(d[4]))
    return Hcs


def get_sqc(apo,achain):
    '''
    Retrieve sequence conservations.
    '''
    # read sequence conservation
    data = open('%s%s.sqc' % (apo,achain))
    D = data.readlines()
    data.close()

    Sqc = {}
    for d in D:
        d = d.strip().split()
        if len(d)>0:
            if d[1]=='X': res, resid = 'UNK', int(d[0])
            else: res, resid = PDB.Polypeptide.one_to_three(d[1]), int(d[0])
            Sqc[(res,resid)] = float(d[2])

    return Sqc



def gather_features(pdb,PDBChainOrder):
    '''
    Gather bioinformatics features (no neighborhood yet).
    '''

    sasa = get_sas(pdb,prog='MOD',probe=3.0)
    prta = get_prt(pdb)
    cnca, cncfa = get_cnc(pdb)
    cvxa = get_cvx(pdb)
    hcsa = {}
    sqca = {}

    for chain in PDBChainOrder:
        hcsa[chain] = get_hcs(pdb,chain)
        sqca[chain] = get_sqc(pdb,chain)
    for chain in PDBChainOrder:
        outf = open(pdb+chain+'.feat','w')
        data = open('%s.pdb' % (pdb))
        D = data.readlines()
        data.close()

        for d in D:
            d = d.strip()
            if d[:4]=='ATOM' and (d[21]==chain):
                atom, res, resid, cid = d[12:16], d[17:20], int(d[22:26]), d[21]
                if atom==' OXT': continue
                sasai = str(sasa[(atom,res,resid,cid)])
                prtai = str(prta[(atom,res,resid,cid)])
                cncai = str(cnca[(res,resid,cid)])
                cvxai = str(cvxa[(atom,res,resid,cid)])
                hcsai = [str(i) for i in hcsa[chain][(res,str(resid))]]
                sqcai = str(sqca[chain][(res,resid)])
                L = [d,sasai,prtai,cncai,cvxai]+hcsai+[sqcai]+list(cncfa)
                outf.write( '\t'.join(L)+'\n' )
        outf.close()
