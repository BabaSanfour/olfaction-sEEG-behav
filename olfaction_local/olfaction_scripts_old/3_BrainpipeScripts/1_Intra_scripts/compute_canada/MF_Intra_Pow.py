
# coding: utf-8

import sys
import numpy as np
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedShuffleSplit as SSS, cross_val_score
from mlxtend.feature_selection import ExhaustiveFeatureSelector as EFS
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from itertools import product

clf = LDA()
cv = SSS(n_splits=5)
cv_inner = SSS(n_splits=4)
feature_names = ('theta', 'alpha', 'beta', 'gamma1', 'gamma2')

PATH = "/home/alsaive/projects/def-kjerbi/alsaive/Intra_EM/"

def MF_Selfeat(su,elec_num,ti):
    pow_data_elec = []
    for power in pow_list:
        dat = power[2:7,elec_num,ti].swapaxes(0,1)
        pow_data_elec.append(dat)
    subsampled_data = np.concatenate(pow_data_elec, axis=0)
    labels = np.hstack([np.array([i]*power.shape[0]) for i, power in enumerate(pow_data_elec)])

    scores, id_sel = [], []
    for inner_index, outer_index in cv.split(subsampled_data, labels):
        inner_data = subsampled_data[inner_index,:]
        inner_labels = labels[inner_index]
        outer_data = subsampled_data[outer_index,:]
        outer_labels = labels[outer_index]

        efs = EFS(clf, 
           min_features=1,
           max_features=5,
           scoring='roc_auc',
           print_progress=False,
           cv=cv_inner,
        )

        efs = efs.fit(inner_data, inner_labels,custom_feature_names=feature_names)
        best_features_index = efs.best_idx_
        best_features_names = efs.best_feature_names_
        clf.fit(inner_data[:, best_features_index], inner_labels)
        scores.append(clf.score(outer_data[:, best_features_index], outer_labels))
        id_sel.append(best_features_names)
    scores = np.asarray(scores)
    id_sel = np.asarray(id_sel)
    np.savez(PATH+'MF_resultsR/'+"{}_{}_{}_Pow_MF_scores.npz".format(su,elec_num,ti), 
        score=scores, id_sel=id_sel,labels=labels)


if __name__ == "__main__":
    
    su = sys.argv[1:][0]

    conds = ['low','high']

    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(PATH+'dataR/'+su+'_odor_'+conds[0]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat0['xpow'][:,:,17:52,:])
    nelecs = mat0['xpow'].shape[1]
    mat1 = np.load(PATH+'dataR/'+su+'_odor_'+conds[1]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat1['xpow'][:,:,17:52,:])

    Parallel(n_jobs=-1)(
        delayed(MF_Selfeat)(su,elec_num,ti) for elec_num,ti in product(range(nelecs),range(35))
        )

