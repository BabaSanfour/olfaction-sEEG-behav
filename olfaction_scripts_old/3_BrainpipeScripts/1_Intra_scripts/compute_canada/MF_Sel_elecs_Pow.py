
# coding: utf-8

import sys
import numpy as np
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedShuffleSplit as SSS, cross_val_score
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from itertools import product
from os import path

clf = LDA()
cv = SSS(n_splits=5)
cv_inner = SSS(n_splits=4)

PATH = "/home/alsaive/projects/def-kjerbi/alsaive/Intra_EM/"

def MF_Selfeat(su,freq,ti):
    savename = "{}/{}_{}_{}_MVPA_scores.npz".format('MF_resultsR',su,ti,freq)

    if path.exists(savename):
        pass
    else:
        pow_data_elec = []
        for power in pow_list:
            dat = power[freq,:,ti].swapaxes(0,1)
            pow_data_elec.append(dat)
        subsampled_data = np.concatenate(pow_data_elec, axis=0)
        labels = np.hstack([np.array([i]*power.shape[0]) for i, power in enumerate(pow_data_elec)])

        scores, id_sel = [], []
        for inner_index, outer_index in cv.split(subsampled_data, labels):
            inner_data = subsampled_data[inner_index,:]
            inner_labels = labels[inner_index]
            outer_data = subsampled_data[outer_index,:]
            outer_labels = labels[outer_index]

            sfs = SFS(
                    clf,
                    k_features=(1, 50),
                    forward=True,
                    floating=False,
                    verbose=0,
                    scoring="roc_auc",
                    cv=cv_inner,
                )

            sfs = sfs.fit(inner_data, inner_labels)
            best_features_index = sfs.k_feature_idx_
            clf.fit(inner_data[:, best_features_index], inner_labels)
            scores.append(clf.score(outer_data[:, best_features_index], outer_labels))
            id_sel.append(best_features_index)
        scores = np.asarray(scores)
        id_sel = np.asarray(id_sel)
        np.savez(PATH+savename, score=scores, id_sel=id_sel, y=labels)


if __name__ == "__main__":
    
    su = sys.argv[1:][0]

    conds = ['low','high']

    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(PATH+'dataR/'+su+'_odor_'+conds[0]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat0['xpow'][:,:,17:52,:])
    nelecs, nfreqs = mat0['xpow'].shape[1], range(2,7)
    mat1 = np.load(PATH+'dataR/'+su+'_odor_'+conds[1]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat1['xpow'][:,:,17:52,:])

    Parallel(n_jobs=-1)(
        delayed(MF_Selfeat)(su,freq,ti) for ti,freq in product(range(35),range(2,7))
        )

