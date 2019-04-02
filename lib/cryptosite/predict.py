#!/usr/bin/env python

"""Do the final prediction of binding site given all features."""

from __future__ import print_function, absolute_import
import pickle
import os
import optparse
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
    import numpy as np
    from sklearn.metrics import confusion_matrix

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

    scaler_pkl = {'linear':'LinearScaler_Final.pkl',
                  'poly':'PolyScaler_Final.pkl',
                  'final':'Scaler_Final_Final.pkl'}[model]
    outmodel_pkl = {'linear':'LinearSVC_FinalModel.pkl',
                    'poly':'PolySVC_FinalModel.pkl',
                    'final':'SVM_Final_Final.pkl'}[model]

    print('Scaling ...')
    with open(os.path.join(cryptosite.config.datadir, scaler_pkl)) as fh:
        scaler = pickle.load(fh)
    X_learn = scaler.transform(X_learn)

    with open(os.path.join(cryptosite.config.datadir, outmodel_pkl)) as fh:
        learner = pickle.load(fh)

    print('Predicting ...')

    # Set _gamma explicitly (earlier versions of cryptosite relied on a hacked
    # local copy of sklearn that did this)
    learner._gamma = 1.0 / X_learn.shape[1]

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

    suffix = {'linear':'lin', 'poly':'pol', 'final':'pol'}[model]
    outn = open(pdb+'.%s.pred' % suffix, 'w')

    print('Writing output files ...')
    outn.write('\t'.join(['PDBID','Res','ResID']+Header+['CryptositeValue'])+'\n')
    for x in range(len(Y_PRED_PROB_ALL)):
        outn.write( '\t'.join(list(NewIndeces[x])+[str(i) for i in X_learn[x]]+[str(Y_PRED_PROB_ALL[x])])+'\n' )
    outn.close()

    write_pdb(pdb, model)

    print('Done!')


def write_pdb(pdb, model='linear'):
    suffix = {'linear':'lin', 'poly':'pol', 'final':'pol'}[model]

    with open(pdb+'.%s.pred' % suffix) as data:
        D = data.readlines()
    out = open('%s.%s.pred.pdb' % (pdb, suffix), 'w')

    Data = {}
    for d in D:
        d = d.strip().split()
        Data[(d[1],d[2])] = ('0.0',d[-1])

    with open('%s_mdl.pdb' % pdb.split('/')[-1]) as data:
        D = data.readlines()

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

def parse_args():
    usage = """%prog [opts] <model_name>

Do the final prediction of binding site given all features.

<model_name> should be the name of the model. The model's 3D structure,
<model_name>_mdl.pdb, and the features file, <model_name>.features,
are read in from the current directory.

Two files are generated on successful prediction:
<model_name>.pol.pred: a simple tab-separated file listing the value of
              the CryptoSite score for each residue.
<model_name>.pol.pred.pdb: a PDB file with the CryptoSite score in the
              occupancy column, for visualization.
"""
    parser = optparse.OptionParser(usage)
    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    model_name = parse_args()
    predict(model_name + '.features', model='final')

if __name__ == '__main__':
    main()
