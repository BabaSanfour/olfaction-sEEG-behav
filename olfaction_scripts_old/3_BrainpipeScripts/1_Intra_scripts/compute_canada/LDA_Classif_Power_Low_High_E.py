# coding: utf-8

from os import path
import sys
from itertools import product
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedKFold as SKFold
from sklearn.metrics import roc_auc_score
from numpy.random import permutation
from joblib import Parallel, delayed

clf = LDA()
cv = SKFold(n_splits=5,shuffle=True)
nperm = 1000

PATH_0 = "/home/alsaive/projects/def-kjerbi/alsaive/Intra_EM/"
PATH = PATH_0+"0_Power_Encoding_EpiPerf_LowHigh/"
SAVE_PATH = PATH_0+"LDA_Power_E_EpiPerf_LowHigh/"

def LDA_E_R_Pow(su,elec,freq):
    elec_label, freq_name = names[elec], freq_names[freq]
    name_auc = (SAVE_PATH+freq_name+'/'+su +'_auc_low_high_'+str(elec_label)+'_('+str(elec)+').npy')
    name_perm = (SAVE_PATH+freq_name+'/'+su +'_perm_'+str(elec_label)+'_('+str(elec)+').npy')
                
    if not path.exists(name_auc):
        pow_data_elec = []
        for i,power in enumerate(pow_list):
            pow_data_elec.append(power[freq,elec].swapaxes(0,1))

        x = np.concatenate(pow_data_elec, axis=0)
        y = np.hstack([np.array([i]*len(power)) for i, power in enumerate(pow_data_elec)])

        auc = np.array([])
        for t in range(x.shape[1]):
            X = x[:,t]
            X = X.reshape(-1, 1)
            score_rep = []
            for i in range(10):
                score_cv = []
                for train_index, test_index in cv.split(X, y):
                    X_train, X_test = X[train_index], X[test_index]
                    y_train, y_test = y[train_index], y[test_index]
                    clf.fit(X=X_train, y=y_train)
                    y_pred = clf.predict(X_test)
                    score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
                score_rep.append(np.mean(score_cv))
            score_rep = np.asarray(score_rep).reshape(1,len(score_rep))
            auc = np.vstack((auc, score_rep)) if np.size(auc) else score_rep
        auc = np.swapaxes(auc,0,1)

        perm_scores = np.array([])
        for t in range(x.shape[1]):
            X = x[:,t]
            X = X.reshape(-1, 1)
            perm_rep = []
            for perm in range(nperm):
                y_perm = y[permutation(len(y))]
                score_cv = []
                for train_index, test_index in cv.split(X, y_perm):
                    X_train, X_test = X[train_index], X[test_index]
                    y_train, y_test = y_perm[train_index], y_perm[test_index]
                    clf.fit(X=X_train, y=y_train)
                    y_pred = clf.predict(X_test)
                    score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
                perm_rep.append(np.mean(score_cv))
            perm_rep = np.asarray(perm_rep).reshape(1,len(perm_rep))
            perm_scores = np.vstack((perm_scores, perm_rep)) if np.size(perm_scores) else perm_rep
        np.save(name_auc, auc)
        np.save(name_perm, perm_scores)


if __name__ == "__main__":
    
    su = sys.argv[1:][0]
    conds = ['low','high']
    freq_names = ['alpha']#['theta','alpha','beta','gamma']
    
    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(PATH+su+'_odor_'+conds[0]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat0['xpow'][:,:,17:52,:])

    #nelecs = mat0['xpow'].shape[1]
    elec = 4
    names = mat0['Mai_RL']

    mat1 = np.load(PATH+su+'_odor_'+conds[1]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat1['xpow'][:,:,17:52,:])
    
    Parallel(n_jobs=-1)(
        delayed(LDA_E_R_Pow)(su,elec,freq) for elec,freq in product(range(nelecs),range(4))
        )