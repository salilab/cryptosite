from __future__ import print_function, absolute_import
import glob
import numpy
from scipy import spatial

def res_packing(pdb, res_list, radius=4.):
    '''
    Calculate packing of the residues. Packing is defined as number of atoms
    within a certain distance from all atoms of a given residues, normalized
    by number of atoms in that residue.
    Parameters:
      - radius: count atoms within this radius
    '''

    data = open('%s.pdb' % pdb)
    D = data.readlines()
    data.close()

    RIDs = dict([(i[1],i[0]) for i in res_list.keys()])

    RRR = {}
    XYZ = []
    for d in D:
        if 'ATOM'==d[:4]:
            if d[17:20]!= 'HOH':
                res, rsid = d[17:20],int(d[22:26])
                x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
                XYZ.append(numpy.array([x,y,z]))
    XYZ = numpy.array(XYZ)
    for r in res_list:
        if r not in RRR: RRR[r]=0
        ZZZ = []
        for a in res_list[r]:
            distances = numpy.sqrt(numpy.sum((XYZ-numpy.array(a))**2,axis=1))
            ZZZ += list(numpy.argwhere(distances<=radius).transpose()[0])
            RRR[r] += len(numpy.argwhere(distances<=radius))
        RRR[r] = len(set(ZZZ))
    Packing = {}
    for r in res_list:
        Packing[r] = RRR[r] / float(len(res_list[r]))
    return Packing

def neighborhood(res_list, radius=4.):
    '''
    Find structural neighbors of a given residue.
    Parameters:
      - radius: maximum distance of a neighbor
    '''

    Graph = {}
    Idx = {}
    A = res_list.keys()

    x=0
    XYZ = []
    for a in A:
        Graph[a] = {}
        for i in res_list[a]:
            Idx[x] = a
            x+=1
            XYZ.append(numpy.array(i))
    XYZ = numpy.array(XYZ)
    Dists = spatial.distance.cdist(XYZ,XYZ)

    C = numpy.argwhere(Dists<=radius)
    for d in C:
        i,j = Idx[d[0]], Idx[d[1]]
        if i!=j: Graph[i][j] = 0

    return Graph



def distance_from_surface(all_atoms, surf_atoms):
    '''
    Calculate distance of a given residue to the surface.
    '''
    #TODO: if I decide to calculate surface using Modeller, then this feature
    # will change as well

    XYZs = []
    Ids = {}
    xs = 0
    for r in surf_atoms:
        for a in surf_atoms[r]:
            XYZs.append(numpy.array(a))
            Ids[xs]=r
            xs+=1
    XYZs = numpy.array(XYZs)

    XYZa = []
    Ida, xa={},0
    for r in all_atoms:
        for a in all_atoms[r]:
            XYZa.append(numpy.array(a))
            Ida[xa]=r
            xa+=1
    XYZa = numpy.array(XYZa)

    Dist = spatial.distance.cdist(XYZs,XYZa)
    Dist2Surf = {}
    for a in range(numpy.shape(Dist)[1]):
        for s in range(numpy.shape(Dist)[0]):
            if Ida[a] not in Dist2Surf: Dist2Surf[Ida[a]] = Dist[s,a]
            else: Dist2Surf[Ida[a]] = min([Dist2Surf[Ida[a]], Dist[s,a]])
    return Dist2Surf



def charge_density(pdb, res_list, radius=5.):
    '''
    Compute charge density in a sphere with a <radius>.
    '''

    charge = {'ALA':0,
            'ARG':1,
            'ASN':0,
            'ASP':-1,
            'CYS':0,
            'GLN':0,
            'GLU':-1,
            'GLY':0,
            'HIS':0.5,
            'ILE':0,
            'LEU':0,
            'LYS':1,
            'MET':0,
            'PHE':0,
            'PRO':0,
            'SER':0,
            'THR':0,
            'TRP':0,
            'TYR':0,
            'VAL':0}


    data = open('%s.pdb' % pdb)
    D = data.readlines()
    data.close()

    XYZ = {}
    for d in D:
        if 'ATOM'==d[:4]:
            if d[17:20]!= 'HOH':
                aid, res, rsid, chainid = d[13:17], d[17:20],int(d[22:26]),d[21]
                x,y,z = float(d[30:38]), float(d[38:46]), float(d[46:54])
                XYZ[(aid, res, rsid, chainid)] = numpy.array([x,y,z])

    data = open('%s.sas' % pdb)
    D = data.readlines()
    data.close()

    SurfRes = {}
    for d in D:
        if 'ATOM'==d[:4]:
            if d[17:20]!= 'HOH':
                aid, res, rsid, chainid = d[13:17], d[17:20],int(d[22:26]),d[21]
                if float(d[60:66]) > 2.: SurfRes[(res,rsid, chainid)] = 0
    for d in D:
        if 'ATOM'==d[:4]:
            if d[17:20]!= 'HOH':
                aid, res, rsid, chainid = d[13:17], d[17:20],int(d[22:26]),d[21]
                if (res,rsid, chainid) in SurfRes: SurfRes[(res,rsid,chainid)] += float(d[60:66])

    CRDS,RS = [],{}
    i = 0
    for xyz in XYZ:
        if (xyz[1],xyz[2],xyz[3]) in SurfRes:
            CRDS.append(XYZ[xyz])
            RS[i] = (xyz[1],xyz[2],xyz[3])
            i+=1
    CRDS = numpy.array(CRDS)

    ChargeDensity = {}
    for r in res_list:
        contacts = {}
        for a in res_list[r]:
            distances = numpy.sqrt(numpy.sum((CRDS-numpy.array(a))**2,axis=1))
            L = numpy.argwhere(distances<=radius).transpose()[0]
            for l in L: contacts[RS[l]] = 0
        C,S=0.,0.
        for i in contacts.keys():
            C += charge[i[0]]
            S += SurfRes[i]
        if S>0.: ChargeDensity[r] = C/S
        else: ChargeDensity[r] = 0.
    return ChargeDensity


def get_patchmap(apo):
    '''
    Retrieve PatchMap feature.
    '''

    # read PatchMap output
    data = open('%s.pdb.ptm' % (apo))
    D = data.readlines()
    data.close()

    Ptm = {}
    for d in D:
        d = d.strip().split()
        if len(d)>0:
            if 'X' in d[1]: res, resid, chainid = 'UNK', int(d[2]), d[3]
            else: res, resid, chainid = d[1], int(d[2]), d[3]
            Ptm[(res,resid,chainid)] = float(d[4])
    return Ptm


# --- APO-HOLO features reader
def res_parser(fil):
    '''
    Parse and gather all BMI features per residue.
    '''

    Header = 'ApoID\tRes\tResID\tChainID\t'
    Header+= '\t'.join(['SAS','PRT','CVX','CNC','SSE','HYD','CHR','SQC','PCK','D2S','PTM']+['NBG','SASn','PRTn','CVXn','CNCn','Un', 'Bn', 'En', 'Gn', 'Hn', 'Sn', 'Tn', 'In','HYDn','CHRn','SQCn','PCKn'])

    # -- output file
    out = open('%s.bmiftr' % (fil,), 'w')
    out.write(Header+'\n')

    # -- read the partial feature list
    apo=fil

    chainfiles = glob.glob(fil[:-1]+'*.feat')
    Apo = {}
    for chainfile in chainfiles:

        chain = chainfile.split('.')[0][-1]
        data = open(chainfile)
        Apo[chain] = data.readlines()
        data.close()

    # -- APO features
    APO = {'AllAtoms':{},
           'SurfAtoms':{},
           'SurfaceArea':{},
           'Protrusion':{},
           'Concavity' : {},
           'Convexity' : {},
           'Rigidity': {},
           'SSE': {'-':0, 'B':0, 'E':0, 'G':0, 'H':0, 'S':0, 'T':0, 'I':0},
           'Hidrophobicity':{},
           'Charge':{},
           'SeqCon':{},
           'PatchMap':{}}

    # - seq. con. normalization factors
    chainfiles = glob.glob(fil[:3]+'*.sqc')
    seqcons = {}
    for chainfile in chainfiles:
        chain = chainfile.split('.')[0][-1]
        data = open(chainfile)
        D = [float(i.strip().split()[-1]) for i in data.readlines()]
        data.close()
        sc_mean, sc_std = numpy.mean(D), numpy.std(D)
        seqcons[chain] = (sc_mean, sc_std)
        D = None

    RESAPO = {}

    # - gather partial features
    for chain in Apo:
        for a in Apo[chain]:
            a = a.strip()
            d = a.split()
            p = (a[17:20], int(a[22:26]), a[21])
            RESAPO[p] = {}
            sas,prt,cnc,cvx,sse,hidro,charge,sqc = d[-10:-2]

            if float(sas) >= 2.:
                if p not in APO['SurfAtoms']: APO['SurfAtoms'][p] = []
                APO['SurfAtoms'][p].append( (float(a[30:38]), float(a[38:46]), float(a[46:54])) )

            if float(sas) >= -10.:
                if p not in APO['AllAtoms']: APO['AllAtoms'][p] = []
                APO['AllAtoms'][p].append( (float(a[30:38]), float(a[38:46]), float(a[46:54])) )

            if p not in APO['SurfaceArea']: APO['SurfaceArea'][p]= [float(sas)]
            else: APO['SurfaceArea'][p].append( float(sas) )

            if p not in APO['Protrusion']: APO['Protrusion'][p] = [float(prt)]
            else: APO['Protrusion'][p].append(float(prt))

            if p not in APO['Convexity']:
                APO['Concavity'][p] = float(cnc)
                APO['Convexity'][p] = float(cvx)
                APO['SSE'][p] = sse
                APO['Hidrophobicity'][p] = float(hidro)
                APO['Charge'][p] = float(charge)
                if seqcons[chain][1]!=0.: APO['SeqCon'][p] = (float(sqc)-seqcons[chain][0])/seqcons[chain][1]
                else: APO['SeqCon'][p] = 0.
                APO['Rigidity'][p] = float(a[60:66])

    # evaluate packing
    Packing = res_packing(apo, APO['AllAtoms'])

    # neighborhood
    Neighbors = neighborhood(APO['AllAtoms'])

    # distance from surface
    Dist2Surface = distance_from_surface(APO['AllAtoms'], APO['SurfAtoms'])

    # evaluate charge density
    ChargeDensity = charge_density(apo, APO['AllAtoms'])

    # get PatchMap feature
    PatchMap = get_patchmap(apo)

    for p in RESAPO.keys():
        RESAPO[p]['SAS'] = numpy.mean(APO['SurfaceArea'][p])
        RESAPO[p]['PRT'] = numpy.mean(APO['Protrusion'][p])
        RESAPO[p]['CVX'] = numpy.mean(APO['Convexity'][p])
        RESAPO[p]['CNC'] = numpy.mean(APO['Concavity'][p])
        RESAPO[p]['SSE'] = APO['SSE'][p] if APO['SSE'][p] !='-' else 'U'
        RESAPO[p]['HYD'] = APO['Hidrophobicity'][p]
        RESAPO[p]['CHR'] = ChargeDensity[p]
        RESAPO[p]['SQC'] = APO['SeqCon'][p]
        RESAPO[p]['PCK'] = Packing[p]
        RESAPO[p]['D2S'] = Dist2Surface[p]
        RESAPO[p]['PTM'] = PatchMap[p]


    Hdr = ['SAS','PRT','CVX','CNC','SSE','HYD','CHR','SQC','PCK','D2S','PTM']
    Hdr+= ['NBG','SASn','PRTn','CVXn','CNCn','Un', 'Bn', 'En', 'Gn', 'Hn', 'Sn', 'Tn', 'In','HYDn','CHRn','SQCn','PCKn']
    # NBG ... number of neighbors

    for p in RESAPO:
        Attn = [len(Neighbors[p])**2] + [0.]*16
        for n in Neighbors[p]:
            if n not in RESAPO: continue
            Attn[1] += RESAPO[n]['SAS']
            Attn[2] += RESAPO[n]['PRT']
            Attn[3] += RESAPO[n]['CVX']
            Attn[4] += RESAPO[n]['CNC']

            if RESAPO[n]['SSE']=='U': Attn[5] += 1
            if RESAPO[n]['SSE']=='B': Attn[6] += 1
            if RESAPO[n]['SSE']=='E': Attn[7] += 1
            if RESAPO[n]['SSE']=='G': Attn[8] += 1
            if RESAPO[n]['SSE']=='H': Attn[9] += 1
            if RESAPO[n]['SSE']=='S': Attn[10] += 1
            if RESAPO[n]['SSE']=='T': Attn[11] += 1
            if RESAPO[n]['SSE']=='I': Attn[12] += 1

            Attn[13] += RESAPO[n]['HYD']
            Attn[14] += RESAPO[n]['CHR']
            Attn[15] += RESAPO[n]['SQC']
            Attn[16] += RESAPO[n]['PCK']

        Attn = numpy.array(Attn)/len(Neighbors[p])

        Att = [apo]+list(p)+[RESAPO[p][i] for i in Hdr[:11]]
        out.write( '\t'.join([str(i) for i in Att+list(Attn)]) + '\n')

    out.close()
