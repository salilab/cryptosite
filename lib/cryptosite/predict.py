#!/usr/bin/env python

"""Do the final prediction of binding site given all features."""

from __future__ import print_function, absolute_import
import numpy as np
from scipy import cluster
import pickle, sys
import os

from sklearn.cluster import KMeans
from sklearn.metrics import confusion_matrix
from sklearn.naive_bayes import GaussianNB
from sklearn import svm
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.decomposition import PCA
from sklearn import preprocessing
from itertools import product
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression
from operator import itemgetter
import cryptosite.config

def get_matrix(inputdata, model='linear'):

    Res = {'CYS': (0, 0, 1, 0, 0), 'ASP': (0, 0, 0, 1, 1), 'SER': (0, 1, 1, 1, 1), 'GLN': (0, 0, 1, 0, 1), 'LYS': (0, 1, 0, 1, 1), 'ILE': (0, 1, 0, 0, 1), 'PRO': (0, 1, 1, 1, 0), 'THR': (1, 0, 0, 0, 0), 'PHE': (0, 1, 1, 0, 1), 'ALA': (0, 0, 0, 0, 0), 'GLY': (0, 0, 1, 1, 1), 'HIS': (0, 1, 0, 0, 0), 'GLU': (0, 0, 1, 1, 0), 'LEU': (0, 1, 0, 1, 0), 'ARG': (0, 0, 0, 0, 1), 'TRP': (1, 0, 0, 0, 1), 'VAL': (1, 0, 0, 1, 1), 'ASN': (0, 0, 0, 1, 0), 'TYR': (1, 0, 0, 1, 0), 'MET': (0, 1, 1, 0, 0)}
    Res = Res.keys()
    SSE = ['B','E','G','H','I','S','T','U']

    data = open(inputdata)
    D = data.readlines()
    data.close()

    Header = D[0].strip().split()

    if model=='poly':
        # the bottom visited for poly SVM
        visited = [Header.index('CNC_mean_300'), Header.index('SQC'), Header.index('D2S')]
        visited+= [Header.index('SQCn'), Header.index('PCKn'), Header.index('Hn')]
        visited+= [Header.index('CN5_std_450'), Header.index('CN5_std_300'), Header.index('CN5_std_350')]
        visited+= [Header.index('CNC'), Header.index('PRT_std_450'), Header.index('CN5_mean_500')]
        visited+= [Header.index('Bn'), Header.index('CHRn'), Header.index('In')]
        visited+= [Header.index('CNC_std_300'), Header.index('CNS_300'),Header.index('SAS14_std_400')]
        visited+= [Header.index('SASn')]

    elif model=='linear':
        # for linear SVM
        visited = [Header.index('CNC_mean_300'), Header.index('SQC'), Header.index('CN5_std_450')]
        visited+= [Header.index('D2S'), Header.index('CNS_300'), Header.index('Hn')]
        visited+= [Header.index('CN5_mean_450'), Header.index('CN5_std_300'), Header.index('SQCn')]
        visited+= [Header.index('CNC_std_350'),Header.index('CNCn'), Header.index('CVX_mean_450')]
        visited+= [Header.index('In')]

    elif model=='final':
        visited = [Header.index('CNC_mean_'), Header.index('SQC'), Header.index('PTM')]

    else:
        print('Wrong model: ', model)
        exit()

    M = []
    Indeces,cnt = {},0
    for d in D[1:]:
        d = d.strip().split('\t')
        LA = []
        for hd in range(len(Header)):
            if hd not in visited: pass
            else:
                if hd==1:
                    L = [0.]*len(Res)
                    L[Res.index(d[1])]=1.
                    LA += L
                elif hd in range(4,8): LA += [float(d[hd])]
                elif hd==8:
                    s = [0.]*len(SSE)
                    s[SSE.index(d[8])]=1.
                    LA += s
                else:
                    LA += [float(d[hd])]
        LA += [float(d[-1])]

        M.append([d[0]]+LA)

        Indeces[cnt] = tuple(d[:3])
        cnt += 1
    return M, [Header[j] for j in sorted(visited)], Indeces



def predict(inputdata, model='linear'):

    print('Reading in the data ...')
    M, Header, Indeces = get_matrix(inputdata, model)

    print('Processing the data for model %s ...' % model.upper())
    pdb=inputdata.split('.')[0]
    print(pdb)
    NewIndeces, newcnt = {},0
    X_learn, Y_learn = [],[]

    for r,m in enumerate(M):
        if len(np.argwhere(np.isnan(np.array(m[1:-1]))==True))>0:
            print(r,m[0])
            print(peter)

        X_learn.append(np.array(m[1:-1]))
        Y_learn.append(m[-1])
        NewIndeces[newcnt] = Indeces[r]
        newcnt+=1

    X_learn = np.array(X_learn)
    X_learn = np.vstack(( X_learn[:,0], X_learn[:,2], X_learn[:,1] )).T

    if model=='linear':

        print('Scaling ...')
        out1 = open(os.path.join(cryptosite.config.datadir,
                                 'LinearScaler_Final.pkl'))
        scaler = pickle.load(out1)
        out1.close()
        X_learn = scaler.transform(X_learn)

        outmodel = open(os.path.join(cryptosite.config.datadir,
                                     'LinearSVC_FinalModel.pkl'))
        learner = pickle.load(outmodel)
        outmodel.close()

    elif model=='poly':

        print('Scaling ...')
        out1 = open(os.path.join(cryptosite.config.datadir,
                                 'PolyScaler_Final.pkl'))
        scaler = pickle.load(out1)
        out1.close()
        X_learn = scaler.transform(X_learn)

        outmodel = open(os.path.join(cryptosite.config.datadir,
                                     'PolySVC_FinalModel.pkl'))
        learner = pickle.load(outmodel)
        outmodel.close()

    elif model=='final':
        print('Scaling ...')
        out1 = open(os.path.join(cryptosite.config.datadir,
                                 'Scaler_Final_Final.pkl'))
        scaler = pickle.load(out1)
        out1.close()
        X_learn = scaler.transform(X_learn)

        outmodel = open(os.path.join(cryptosite.config.datadir,
                                     'SVM_Final_Final.pkl'))
        learner = pickle.load(outmodel)
        outmodel.close()

    else:
        print('Unknown model: ', model)
        print(peter)


    print('Predicting ...')
    Y_pred = learner.predict(X_learn)
    CM = confusion_matrix(Y_learn,Y_pred)
    print()
    print("Confusion matrix for: ",pdb)
    print(CM)
    print()

    # output
    Y_pred_prob = learner.predict_proba(X_learn)
    Y_PRED_PROB_ALL = list(Y_pred_prob[:, 1])
    Y_TEST_ALL = list(Y_learn)

    if model=='linear': outn = open(pdb+'.lin.pred','w')
    elif model=='poly' or model=='final': outn = open(pdb+'.pol.pred','w')
    else: print(peter)

    print('Writing output files ...')
    outn.write('\t'.join(['PDBID','Res','ResID']+Header+['CryptositeValue'])+'\n')
    for x in range(len(Y_PRED_PROB_ALL)):
        outn.write( '\t'.join(list(NewIndeces[x])+[str(i) for i in X_learn[x]]+[str(Y_PRED_PROB_ALL[x])])+'\n' )
    outn.close()

    if model=='linear': write_pdb(pdb,model='linear')
    elif model=='poly' or model=='final': write_pdb(pdb,model='poly')
    else: print(peter)

    print('Done!')


def write_pdb(pdb,model='linear'):

    if model=='linear':
        data = open(pdb+'.lin.pred')
        D = data.readlines()
        data.close()
        out = open('%s.lin.pred.pdb' % pdb, 'w')
    elif model=='poly':
        data = open(pdb+'.pol.pred')
        D = data.readlines()
        data.close()
        out = open('%s.pol.pred.pdb' % pdb, 'w')
    else: print(peter)

    Data = {}
    for d in D:
        d = d.strip().split()
        Data[(d[1],d[2])] = ('0.0',d[-1])

    data = open('%s_mdl.pdb' % pdb.split('/')[-1])
    D = data.readlines()
    data.close()

    for d in D:
        if 'ATOM'==d[:4]:
            p = (d[17:20], str(int(d[22:26])))

            try: pred = Data[p]
            except KeyError: pred = ('0.0', '0.0')

            v = '%.2f' % (float(pred[1])*100)
            v = (6-len(v))*' '+v
            line = d[:56]+pred[0]+'0'+v+'\n'
            out.write(line)

        else: out.write(d)
    out.close()

def main():
    predict(sys.argv[-1]+'.features', model='final')

if __name__ == '__main__':
    main()
