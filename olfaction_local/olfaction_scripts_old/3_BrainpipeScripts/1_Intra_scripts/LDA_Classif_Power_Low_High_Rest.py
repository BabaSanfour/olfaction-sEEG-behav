# coding: utf-8

from os import path
from itertools import product
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedKFold as SKFold
from sklearn.metrics import roc_auc_score
from numpy.random import permutation

clf = LDA()
cv = SKFold(n_splits=5,shuffle=True)
nperm = 1000


PATH = "/home/karim/mp2/Intra_EM/0_Power_R_pre_stim_1s/"

def LDA_Rest_Pow(su,elec,freq):
    
    elec_label, freq_name = names[elec], freq_names[freq]
    name_auc = (save_path+freq_name+'/'+su +'_auc_low_high_'+str(elec_label)+'_('+str(elec)+').npy')
    name_perm = (save_path+freq_name+'/'+su +'_perm_'+str(elec_label)+'_('+str(elec)+').npy')
                
    if not path.exists(name_auc):
        pow_data_elec = []
        for i,power in enumerate(pow_list):
            pow_data_elec.append(power[freq,elec][np.newaxis].swapaxes(0,1))

        x = np.concatenate(pow_data_elec, axis=0)
        X = x.reshape(-1, 1)
        y = np.hstack([np.array([i]*len(power)) for i, power in enumerate(pow_data_elec)])

        score_rep = np.array([])
        for i in range(10):
            score_cv = []
            for train_index, test_index in cv.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                clf.fit(X=X_train, y=y_train)
                y_pred = clf.predict(X_test)
                score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
            score_rep = np.vstack((score_rep,np.mean(score_cv)) if np.size(score_rep) else np.mean(score_cv)
        
        perm_rep = np.array([])
        for perm in range(nperm):
            y_perm = y[permutation(len(y))]
            score_cv = []
            for train_index, test_index in cv.split(X, y_perm):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y_perm[train_index], y_perm[test_index]
                clf.fit(X=X_train, y=y_train)
                y_pred = clf.predict(X_test)
                score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
            perm_rep = np.vstack((perm_rep,np.mean(score_cv)) if np.size(perm_rep) else np.mean(score_cv)
        np.save(name_auc, auc)
        np.save(name_perm, perm_scores)
        del X, auc, pow_data_elec


if __name__ == "__main__":
    
    su = sys.argv[1:][0]
    conds = ['low','high']
    freq_names = ['VLFC','delta','theta','alpha','beta','low_gamma','high_gamma']
    
    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(PATH+su+'_odor_'+conds[0]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat0['xpow'])

    nelecs = mat0['xpow'].shape[1]
    names = mat0['Mai_RL']

    mat1 = np.load(PATH+su+'_odor_'+conds[1]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat1['xpow'])
    
    Parallel(n_jobs=-1)(
        delayed(LDA_Rest_Pow)(su,elec,freq) for elec,freq in product(range(nelecs),range(2,7))
        )

