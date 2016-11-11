#!/usr/bin/env python

"""Initial setup and preparation of AllosMod inputs (short version)"""

from __future__ import print_function, absolute_import
import sys
from cryptosite.cleaning import *
from cryptosite.seq_conservation import *
from cryptosite.hyd_chr_sse import *
from cryptosite.bmi_feature_parser import *
from cryptosite.res_parser_bmi import *
from cryptosite.patch_mapper import *
import os
import glob

def main():
    pdb,chains = 'XXX', ''
    paramdata = open('param.txt')
    params = paramdata.readlines()
    paramdata.close()
    chains = params[2].strip().strip(' ').split(',')

    # --- get input PDB,input Seq
    if 'input.pdb' in os.listdir('.'): os.system('cp input.pdb %s.pdb' % (pdb, ))
    else: os.system('cp *.ent %s.pdb' % (pdb, ))



    # --- trim if NMR


    # --- delete unused chains from pdb + extract PDBChainOrder

    PDBChainOrder=[]
    myfile = open('%s.pdb' % pdb,'r')
    data = myfile.read().split('ENDMDL')[0] # trim mulitple models
    file=data.split('\n')
    myfile.close()




    myfile=open('%s.pdb' % pdb,'w')
    for line in file:
        if line[:5]=='ATOM ' or line[:3]=='TER':
            if len(line)<60: continue
            if line[21] in chains:
                if len(line)>77 and (line[77] =='D' or line[77]=='H'): continue
                line=line.strip('/n')
                myfile.write("%s\n" % line)
                if line[21] not in PDBChainOrder:
                    PDBChainOrder.append(line[21])

    myfile.close()


    # --- write seqeunces
    output = open('input.seq','w')
    querySeqs = {}
    for chain in chains:
        querySeq = get_pdb_seq(pdb+'.pdb', chain)
        querySeqs[chain] = querySeq
        output.write('>'+chain+'\n'+querySeq+'\n')
    output.close()



    print('PDB: ',pdb,'\tCHAIN: ',chains, PDBChainOrder)
    # --- start of for loop going through all chains for sequence features

    ChainLenghts, SubjectSeqList, QuerySeqList = {},{},{}

    for chain in PDBChainOrder:
        try: querySeq = querySeqs[chain]
        except KeyError:
            print('CHAIN: ',chain, 'not in the PDB file')
            raise
        print(querySeq, querySeqs[chain])
        ChainLenghts[chain]=(querySeq[0],len(querySeq),querySeq[-1])
        sbjctSeq = ''

        sbjctSeq = querySeq

        # --- refine the structure (Modeller)
        strcsq, seqsq = muscleAlign(querySeq, sbjctSeq, pdb, chain)
        SubjectSeqList[chain] = seqsq
        QuerySeqList[chain] = strcsq

        # --- extract bioinformatics features

        # -- sequence conservation

        output = open('test.seq','w')
        output.write('>%s%sq\n' % (pdb, chain) + sbjctSeq)
        output.close()
        run_blast(pdb+chain)

        parse_blast(pdb+chain+'.blast', pdb+chain, sbjctSeq)



    # --- change ali to pir
    out = open('alignment.pir','w')

    strcsq = '/'.join([QuerySeqList[chain] for chain in PDBChainOrder])
    seqsq ='/'.join([SubjectSeqList[chain] for chain in PDBChainOrder])
    out.write('>P1; %s\nstructure:%s:FIRST    :%s:LAST  :%s:: : :\n' % (pdb,pdb,PDBChainOrder[0],PDBChainOrder[-1]))
    out.write(strcsq+'*\n')

    out.write('\n>P1; %s_X\nsequence:%s_X:    :%s:  :%s:: : :\n' % (pdb.lower(),pdb.lower(),PDBChainOrder[0],PDBChainOrder[-1]))
    out.write(seqsq+'*\n')

    out.close()


    # --- build model of all the chains

    #RESMAP=build_model(pdb, PDBChainOrder, ChainLenghts)

    print("Printing MODELLER input:")
    print(pdb)
    print(PDBChainOrder)
    print(ChainLenghts)

    build_model(pdb, PDBChainOrder)

    # --- Map SeqConservation to new residue numbering calculate HydChrSSE

    mchains=map(chr, range(65, 65+len(PDBChainOrder)))

    L=0
    for chain in PDBChainOrder:
        seqdat=open(pdb+chain+'.sqc')
        sqdat=seqdat.readlines()
        seqdat.close()

        seqdat=open(pdb+'_mdl'+mchains[PDBChainOrder.index(chain)]+'.sqc','w')
        for sql in sqdat:
            sql=sql.split()
            sql[0]=str(int(sql[0])+L)
            seqdat.write('\t'.join(sql)+'\n')

        L+=len(SubjectSeqList[chain])
        seqdat.close()
        # -- hydrophobicity, charge, SSEs
        HydChrSSE(pdb+'_mdl', mchains[PDBChainOrder.index(chain)])


    # --- calculate PatchMap feature
    patchmap_feature(pdb+'_mdl')


    # --- gather residue-based BMI features
    gather_features(pdb+'_mdl',mchains)


    #for chain in PDBChainOrder:
    res_parser(pdb+'_mdl')


    # --- prepare AllosMod file

    f = pdb+'_mdl.bmiftr'
    data = open(f)
    F1 = len(data.readlines())
    data.close()

    data = open('alignment.pir')
    S = data.read()
    D = S.split(':')[-1].replace('\n','')[:-1].replace('-','')
    P = '>'+S.split('>')[1].replace(pdb,pdb+'.pdb').replace('structure','structureX')
    data.close()

    os.system('mkdir -p %s' % pdb)
    out = open('%s/align.ali' % pdb, 'w')
    out.write( P )
    out.write(">P1;pm.pdb\n" )
    out.write("structureX:pm.pdb:1    :%s:%i  :%s::::\n" % (chain,len(D),chain))
    out.write(D+'*\n')
    out.close()

    os.system('cp %s.pdb %s' % (pdb,pdb))

    out1 = open('%s/input.dat' % pdb, 'w')
    out1.write('rAS=1000\nNRUNS=25\nSCRAPP=True\nMDTemp=SCAN')
    out1.close()

    out2 = open('%s/list' % pdb, 'w')
    out2.write('%s.pdb' % pdb)
    out2.close()

if __name__ == '__main__':
    main()
