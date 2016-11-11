from __future__ import print_function, absolute_import
from Bio import PDB
from numpy import linalg, array
import subprocess
import cryptosite.am_bmi

def get_cnc(apo):
    '''
    Find pockets using Fpocket algorithm.
    '''

    subprocess.check_call(["fpocket", "-f", apo+'.pdb'])

    with open('%s_out/%s_out.pdb' % (apo,apo)) as data:
        D = data.readlines()

    Pockets = {}
    for d in D:
        if d[17:20]=='STP' and d[:6]=='HETATM':
            pn = int(d[22:26])
            coords = (float(d[30:38]), float(d[38:46]), float(d[46:54]))
            if pn not in Pockets:
                Pockets[pn] = [coords]
            else:
                Pockets[pn].append(coords)

    PocketInfo = {}
    with open('%s_out/%s_info.txt' % (apo,apo)) as data:
        Z = data.read().split('\n\n')

    for p in Z:
        p = p.split('\n')
        pn = ''
        for i in p:
            if i[:6]=='Pocket':
                pn = int(i.split(':')[0][6:])
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
                    if dist<mini:
                        mini=dist
                if mini<=5.:
                    adist[p] = mini
            atom, res, resid, cid  = d[12:16], d[17:20], int(d[22:26]), d[21]
            Residues[(res,resid,cid)] = 0.

            if len(adist)>0:
                if max([PocketInfo[i]
                       for i in adist])>Residues[(res,resid,cid)]:
                    Residues[(res,resid,cid)] = max([PocketInfo[i]
                                                     for i in adist])
    print(Residues)
    return Residues, ('1','1')

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

    sasa = cryptosite.am_bmi.get_sas(pdb + '.pdb', probe=3.0)
    prta = cryptosite.am_bmi.get_prt(pdb + '.pdb')
    cnca, cncfa = get_cnc(pdb)
    cvxa = cryptosite.am_bmi.get_cvx(pdb + '.pdb')
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
