#!/usr/bin/env python

"""Initial setup and preparation of AllosMod inputs"""

from __future__ import print_function, absolute_import
import sys
import cryptosite.cleaning
import cryptosite.seq_conservation
import cryptosite.hyd_chr_sse
import cryptosite.bmi_feature_parser
import cryptosite.res_parser_bmi
import cryptosite.patch_mapper
import os
import glob
import shutil
import optparse

class MissingChainsError(Exception):
    pass

def extract_chains(fname, chains):
    """Extract the named chains from the file.
       Return a list of chains in the order they appear in the file,
       and modify the file in place to contain only these chains."""
    PDBChainOrder=[]
    with open(fname,'r') as myfile:
        data = myfile.read().split('ENDMDL')[0] # trim mulitple models
        file=data.split('\n')

    with open(fname, 'w') as myfile:
        for line in file:
            if line[:5]=='ATOM ' or line[:3]=='TER':
                if len(line)<60: continue
                if line[21] in chains:
                    if len(line)>77 and (line[77] =='D' or line[77]=='H'):
                        continue
                    line=line.strip('/n')
                    myfile.write("%s\n" % line)
                    if line[21] not in PDBChainOrder:
                        PDBChainOrder.append(line[21])
    missing_chains = set(chains) - set(PDBChainOrder)
    if missing_chains:
        raise MissingChainsError("The following chains were not found in the "
                              "input PDB file: %s" % ", ".join(missing_chains))
    return PDBChainOrder

def setup(fname, chains, short):
    pdb = 'XXX'

    # --- get input PDB,input Seq
    shutil.copy(fname, '%s.pdb' % pdb)

    # --- trim if NMR


    # --- delete unused chains from pdb + extract PDBChainOrder
    PDBChainOrder = extract_chains('%s.pdb' % pdb, chains)

    # --- write sequences
    querySeqs = {}
    with open('input.seq', 'w') as output:
        for chain in chains:
            querySeq = cryptosite.cleaning.get_pdb_seq(pdb+'.pdb', chain)
            querySeqs[chain] = querySeq
            output.write('>'+chain+'\n'+querySeq+'\n')

    print('PDB: ',pdb,'\tCHAIN: ',chains, PDBChainOrder)
    # --- start of for loop going through all chains for sequence features

    ChainLengths, SubjectSeqList, QuerySeqList = {},{},{}

    for chain in PDBChainOrder:
        try: querySeq = querySeqs[chain]
        except KeyError:
            print('CHAIN: ',chain, 'not in the PDB file')
            raise
        print(querySeq, querySeqs[chain])
        ChainLengths[chain]=(querySeq[0],len(querySeq),querySeq[-1])
        sbjctSeq = ''

        sbjctSeq = querySeq

        # --- refine the structure (Modeller)
        strcsq, seqsq = cryptosite.cleaning.muscleAlign(querySeq, sbjctSeq,
                                                        pdb, chain)
        SubjectSeqList[chain] = seqsq
        QuerySeqList[chain] = strcsq

        # --- extract bioinformatics features

        # -- sequence conservation

        with open('test.seq','w') as output:
            output.write('>%s%sq\n' % (pdb, chain) + sbjctSeq)
        cryptosite.seq_conservation.run_blast(pdb+chain)

        cryptosite.seq_conservation.parse_blast(pdb+chain+'.blast', pdb+chain,
                                                sbjctSeq)

    # --- change ali to pir
    with open('alignment.pir', 'w') as out:
        strcsq = '/'.join([QuerySeqList[chain] for chain in PDBChainOrder])
        seqsq ='/'.join([SubjectSeqList[chain] for chain in PDBChainOrder])
        out.write('>P1; %s\nstructure:%s:FIRST    :%s:LAST  :%s:: : :\n'
                  % (pdb, pdb, PDBChainOrder[0], PDBChainOrder[-1]))
        out.write(strcsq+'*\n')

        out.write('\n>P1; %s_X\nsequence:%s_X:    :%s:  :%s:: : :\n'
                  % (pdb.lower(), pdb.lower(), PDBChainOrder[0],
                     PDBChainOrder[-1]))
        out.write(seqsq + '*\n')

    # --- build model of all the chains

    #RESMAP=build_model(pdb, PDBChainOrder, ChainLengths)

    print("Printing MODELLER input:")
    print(pdb)
    print(PDBChainOrder)
    print(ChainLengths)

    cryptosite.cleaning.build_model(pdb, PDBChainOrder)

    # --- Map SeqConservation to new residue numbering calculate HydChrSSE

    mchains=[chr(x) for x in range(65, 65+len(PDBChainOrder))]

    L=0
    for chain in PDBChainOrder:
        with open(pdb + chain + '.sqc') as seqdat:
            sqdat=seqdat.readlines()

        with open(pdb + '_mdl' + mchains[PDBChainOrder.index(chain)]
                  + '.sqc', 'w') as seqdat:
            for sql in sqdat:
                sql=sql.split()
                sql[0]=str(int(sql[0])+L)
                seqdat.write('\t'.join(sql)+'\n')

            L+=len(SubjectSeqList[chain])

        # -- hydrophobicity, charge, SSEs
        cryptosite.hyd_chr_sse.HydChrSSE(pdb+'_mdl',
                                         mchains[PDBChainOrder.index(chain)])


    # --- calculate PatchMap feature
    cryptosite.patch_mapper.patchmap_feature(pdb+'_mdl')


    # --- gather residue-based BMI features
    cryptosite.bmi_feature_parser.gather_features(pdb+'_mdl',mchains)


    #for chain in PDBChainOrder:
    cryptosite.res_parser_bmi.res_parser(pdb+'_mdl')


    # --- prepare AllosMod file

    f = pdb+'_mdl.bmiftr'
    with open(f) as data:
        F1 = len(data.readlines())

    with open('alignment.pir') as data:
        S = data.read()
    D = S.split(':')[-1].replace('\n','')[:-1].replace('-','')
    P = '>'+S.split('>')[1].replace(pdb,pdb+'.pdb') \
                           .replace('structure','structureX')

    if not os.path.exists(pdb):
        os.mkdir(pdb)

    with open('%s/align.ali' % pdb, 'w') as out:
        out.write(P)
        out.write(">P1;pm.pdb\n" )
        out.write("structureX:pm.pdb:1    :%s:%i  :%s::::\n"
                  % (chain,len(D),chain))
        out.write(D+'*\n')

    shutil.copy('%s.pdb' % pdb, pdb)

    with open('%s/input.dat' % pdb, 'w') as out1:
        if short:
            out1.write('rAS=1000\nNRUNS=25\nMDTemp=SCAN\n')
        else:
            out1.write('rAS=1000\nNRUNS=50\nMDTemp=SCAN\n')

    with open('%s/list' % pdb, 'w') as out2:
        out2.write('%s.pdb' % pdb)

def parse_args():
    usage = """%prog [opts] <pdb> <chains>

Do initial setup and preparation of AllosMod inputs.
<pdb> should be the full path to a structure file in PDB format.
<chains> should be a comma-separated list of chain IDs.
"""
    parser = optparse.OptionParser(usage)
    parser.add_option("--short", action="store_true",
                      help="run a quicker AllosMod simulation")
    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    chains = args[1].strip().strip(' ').split(',')
    return args[0], chains, opts.short

def main():
    pdb, chains, short = parse_args()
    setup(pdb, chains, short)

if __name__ == '__main__':
    main()
