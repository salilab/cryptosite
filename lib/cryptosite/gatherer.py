#!/usr/bin/env python

"""Gather all feature information into a single file."""

import sys, glob
import numpy as np
import os
import random

def process_directory(tts):
    pdbs = os.listdir(tts + '/pred_dECALCrAS1000')
    # pdbs = [i for i in pdbs if '_' in i and 'pdb' in i and PDB in i]
    pdbs = [i for i in pdbs if '_' in i and 'pdb' in i ]
    pdbs = list(set(pdbs))

    QM,QD = [],[]
    PM,PD,PM5,PD5,PMs,PDs,Len =[],[],[],[],[],[],[]
    AM = {'sas14':{},'sas30':{},'prt':{},'cvx':{}}
    RES = []

    for pdb in pdbs:
        efiles = glob.glob(tts + '/pred_dECALCrAS1000/'+pdb+'/energy.dat')

        print tts + '/pred_dECALCrAS1000/'+pdb+'/energy.dat'

        for e in efiles:
            pdbf = e.rsplit('/',2)[1]
            #pdb = pdbf.split('_')[0]

            # --- read energy file and determine the temperature
            data = open(e)
            D = data.readlines()
            data.close()
            E = [float(i.strip().split()[-1]) for i in D]
            # --- read in Qi
            try:
                #data = open(e.rsplit('/',1)[0]+'/qioft_pm_'+pdb.split('.')[0]+'.pdb_11sc.dat')
                data = open(e.rsplit('/',1)[0]+'/qioft_pm_XXX.pdb_11sc.dat')
                D = data.readlines()
                data.close()
                q = []
                #for d in D: Q[temp].append( np.array([float(i) for i in d.strip().split()]) )
                if len(D)>0:# and 'THEORETICAL' not in D[0]:
                    for d in D:
                        q.append( np.array([float(i) for i in d.strip().split()]) )
                    q = np.array(q)
                    QM.append( np.mean(q,axis=0) )
                    QD.append( np.std(q,axis=0) )
            except (IOError, ValueError): pass

            # --- read in pockets, --- chainID d[2]?
            try:
                pocket_files = glob.glob(tts + '/pred_dECALCrAS1000/'+pdb+'/pockets.out')
                data = open(pocket_files[0])
                D = data.readlines()
                data.close()
                p = []
                Res = []
                len_limit = len(D[0].strip().split())-2
                for d in D[1:]:
                    d = d.strip().split()
                    if len(d)>5:
                        sd = [float(i) for i in d[3:]]
                        if len(sd)==len_limit: p.append(np.array(sd))
                        elif len(sd)<len_limit: p.append(np.array(sd+[0.]*(len_limit-len(sd))))
                    Res.append((d[0],int(d[1]),d[2]))
                p = np.array(p)
                PM.append( np.mean(p,axis=1) )
                PD.append( np.std(p,axis=1) )
                PMs.append( (p>=0.4).sum(axis=1) / float(np.shape(p)[1]) )
                percentile = int(round(len(Res)*0.05,0))
                PM5.append( np.sort(p,axis=1)[:,-percentile:].mean(axis=1) )
                PD5.append( np.sort(p,axis=1)[:,-percentile:].std(axis=1) )
                if len(RES)==0: RES = Res
            except (IOError, ValueError, IndexError):
                print 'Error'
                print d
                pass
            # --- read in geometry
            files = glob.glob(tts+'/pred_dECALCrAS1000/'+pdb+'/am_features.out')
            if len(files)!=1: continue
            for fil in files:
                data = open(fil)
                D = data.readlines()
                data.close()

                #print pdbf,fil
                for d in D:
                    d = d.strip().split()
                    if len(d)!=7: continue
                    try: res = (d[0],int(d[1]),d[2])
                    except (IndexError, ValueError): continue
                    if res not in AM['sas14']:
                        try:
                            AM['sas14'][res] = [float(d[3])]
                            AM['sas30'][res] = [float(d[4])]
                            AM['prt'][res]   = [float(d[5])]
                            AM['cvx'][res]   = [float(d[6])]
                        except (IndexError, ValueError): pass
                    else:
                        try:
                            AM['sas14'][res].append(float(d[3]))
                            AM['sas30'][res].append(float(d[4]))
                            AM['prt'][res].append(float(d[5]))
                            AM['cvx'][res].append(float(d[6]))
                        except (ValueError, IndexError): pass
            #print pdbf,pdb,temp,percentile,len(files)
    qm = np.array(QM)
    qd = np.array(QD)
    pm = np.array(PM)
    pd = np.array(PD)
    pm5= np.array(PM5)
    pd5= np.array(PD5)
    pms= np.array(PMs)

    qm = np.mean(qm, axis=0)
    qd = np.mean(qd, axis=0)
    pm = np.mean(pm, axis=0)
    pd = np.mean(pd, axis=0)
    pm5= np.mean(pm5,axis=0)
    pd5= np.mean(pd5,axis=0)
    pms= np.mean(pms,axis=0)

    sas14m, sas14d = [],[]
    sas30m, sas30d = [],[]
    prtm, prtd = [],[]
    cvxm, cvxd = [],[]

    out = open(pdb.split('.')[0]+'.am','w')
    H = ['Res','ResID','SAS14_mean_','SAS14_std_','SAS30_mean_','SAS30_std_','PRT_mean_','PRT_std_']
    H += ['CVX_mean_','CVX_std_','QI_mean_','QI_std_','CNC_mean_','CNC_std_']
    H += ['CN5_mean_','CN5_std_','CNS_']
    out.write('\t'.join(H)+'\n')
    for res in RES:
        L = [ res[0],res[1],res[2],np.mean(AM['sas14'][res]),np.std(AM['sas14'][res])]
        L +=[np.mean(AM['sas30'][res]),np.std(AM['sas30'][res])]
        L +=[np.mean(AM['prt'][res]),np.std(AM['prt'][res])]
        L +=[np.mean(AM['cvx'][res]),np.std(AM['cvx'][res])]
        L +=[qm[res[1]-1], qd[res[1]-1], pm[res[1]-1], pd[res[1]-1]]
        L +=[pm5[res[1]-1], pd5[res[1]-1], pms[res[1]-1]]
        out.write('\t'.join([str(i) for i in L])+'\n')
    out.close()
    print np.shape(qm), np.shape(qd), np.shape(pm), np.shape(pd),np.shape(pm5),np.shape(pms) ,pdb


    ### --- create the final file with all the features

    #out = open('1LAYA2.features','w')
    Header = ['ApoID','Res','ResID','ChainID','SAS','PRT','CVX','CNC','SSE','HYD','CHR','SQC','PCK','D2S','PTM','NBG','SASn','PRTn','CVXn','CNCn','Un','Bn','En','Gn','Hn','Sn','Tn','In','HYDn','CHRn','SQCn','PCKn','SAS14_mean_','SAS14_std_','SAS30_mean_','SAS30_std_','PRT_mean_','PRT_std_','CVX_mean_','CVX_std_','QI_mean_','QI_std_','CNC_mean_','CNC_std_','CN5_mean_','CN5_std_','CNS_','CBS']
    #out.write('\t'.join(Header)+'\n')


    files = glob.glob('*.am')
    pdbs = set([i.split('_')[0].split('/')[-1] for i in files])
    for a in pdbs: #if 1:

            #a = '1LAYA2'
        fff = a.split('.')[0]
        data = open(fff+'_mdl.bmiftr')
        D = data.readlines()
        data.close()

        out = open('%s.features' % (fff,),'w')
        out.write('\t'.join(Header)+'\n')

        H = D[0].strip().split()
        A = {}
        for d in D[1:]:
            d = d.strip().split()
            A[(d[1],int(d[2]),d[3])] = d


        data = open('%s.am' % fff)
        D = data.readlines()
        data.close()

        H += ['SAS14_mean_','SAS14_std_','SAS30_mean_','SAS30_std_','PRT_mean_','PRT_std_']
        H += ['CVX_mean_','CVX_std_','QI_mean_','QI_std_','CNC_mean_','CNC_std_']
        H += ['CN5_mean_','CN5_std_','CNS_']
        H += ['CBS']

        for d in D[1:]:
            print d
            d = d.strip().split()
            try:
                A[(d[0],int(d[1]),d[2])] += d[3:]
            except:
                pass


             #out.write('\t'.join(Header)+ '\t'.join(H)+'\n')

        for i in A:
            out.write('\t'.join(A[i]+['0'])+'\n')
        out.close()

    print H, len(H)

def main():
    process_directory('/scrapp/AM/'+sys.argv[-2]) #random filename generated by backend

if __name__ == '__main__':
    main()
