
# coding: utf-8

# ## Import Libraries

# In[2]:


from os import path
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import scipy.io as sio

from brainpipe.classification import *
from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt
from brainpipe.visual import *
from brainpipe.statistics import *
from scipy.stats import *


# ## Power Decoding - Partial//Detailed Encoding
# ### For ALL time points

# In[6]:


PATH ='/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/database/Retrieval_Rest_LowHigh/'
mat = np.load(PATH+'VACJ_odor_high_bipo_sel_physFT.npz')['x']
print(mat.shape)


# In[ ]:


from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedKFold as SKFold
from sklearn.metrics import roc_auc_score
from numpy.random import permutation

conds,phases, subjects = ['low','high'],['odor'],['SEMC','VACJ','FERJ','LEFC','PIRJ','CHAF']#,'MICP','VACJ','SEMC','LEFC','PIRJ']
color_codes = ['darkorange','blue']
st = study('Olfacto')
freqs = 7
nperm = 1000

for su,phase in product(subjects,phases):
    path_pow = path.join(st.path, 'feature/TPSim_Enc_Ret_By_Odor/TPS_R_p_by_cond/')
    save_path = path.join(st.path, 'classified/Classif_TPSim_Encoding_Retrieval_Odor_1000perm/')

    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(path.join(path_pow, su+'_'+phase+'_'+conds[0]+'_bipo_sel_physFT_pow'+bsl+'.npz'))
    names, channels, freq_names, time = mat0['Mai_RL'], mat0['channels'],mat0['fname'], mat0['time']
    #print(mat0['xpow'].shape,time.shape,time[17:42])
    pow_list.append(mat0['xpow'][:,:,18,:][:,:,np.newaxis,:])
    nelecs = mat0['xpow'].shape[1]
    mat1 = np.load(path.join(path_pow, su+'_'+phase+'_'+conds[1]+'_bipo_sel_physFT_pow'+bsl+'.npz'))
    pow_list.append(mat1['xpow'][:,:,18,:][:,:,np.newaxis,:])
    print (su, 'power shape: ', [pow.shape for pow in pow_list])
    # =========================== Select Power for 1 elec 1 freq =================================                 
    iterator = range(nelecs)
    for elec_num in iterator:
        for freq in range(2,freqs):
            elec, elec_label, freq_name = channels[elec_num], names[elec_num], freq_names[freq]
            print ('elec ', elec, 'elec_label ', elec_label)
            #Filenames to save
            name_auc = (save_path+str(freq)+'_'+freq_name+'/auc/'+su +'_auc_'+conds[0]+'_'+conds[1]+'_'+str(elec_label)+'_('+str(elec_num)+').npy')
            name_perm = (save_path+str(freq)+'_'+freq_name+'/auc/'+su +'_perm_'+str(elec_label)+'_('+str(elec_num)+').npy')
                        
            if path.exists(name_auc):
                print(su,bsl,phase,elec_num,freq,'already computed')
            else:
                print('--» processing',su, 'elec', elec_num,'/',nelecs, 'freq',freq)
                pow_data_elec = []
                for i,power in enumerate(pow_list):
                    pow_data_elec.append(power[freq,elec_num].swapaxes(0,1))
                nwin = power.shape[1]

        # =============================  Classification Computation ============================================================           
                # create a data matrix, concatenate along the trial dimension
                x = np.concatenate(pow_data_elec, axis=0)
                print ('Size of the concatenated data: ', x.shape, 'Number time windows : ', x.shape[1])
                y = np.hstack([np.array([i]*len(power)) for i, power in enumerate(pow_data_elec)])
                print ('Size of label for classif: ', len(y))

                auc = np.array([])
                for t in range(x.shape[1]):
                    X = x[:,t]
                    X = X.reshape(-1, 1)
                    score_rep = []
                    for i in range(10):
                        k = 5
                        skf = SKFold(n_splits=k, random_state=None, shuffle=True)
                        skf.get_n_splits(X, y)
                        score_cv = []
                        for train_index, test_index in skf.split(X, y):
                            clf = LDA()
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
                        for train_index, test_index in skf.split(X, y_perm):
                            clf = LDA()
                            X_train, X_test = X[train_index], X[test_index]
                            y_train, y_test = y_perm[train_index], y_perm[test_index]
                            clf.fit(X=X_train, y=y_train)
                            y_pred = clf.predict(X_test)
                            score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
                        perm_rep.append(np.mean(score_cv))
                    perm_rep = np.asarray(perm_rep).reshape(1,len(perm_rep))
                    perm_scores = np.vstack((perm_scores, perm_rep)) if np.size(perm_scores) else perm_rep
                perm_scores = np.swapaxes(perm_scores,0,1)           
                th_0_05_perm = perm_pvalue2level(perm_scores, p=0.01, maxst=True)
                th_0_01_perm = perm_pvalue2level(perm_scores, p=0.001, maxst=True)
                print('th_perm 005: ', th_0_05_perm[0], '001',th_0_01_perm[0], 'auc_max', np.max(auc))

        # ============================== PLOT POWER ANALYSIS + STATS & DECODING ACCURACY ===================================================
#                 # plot and figure parameters
#                 xfmt = ScalarFormatter(useMathText=True)
#                 xfmt.set_powerlimits((0,3))
#                 fig = plt.figure(1,figsize=(7,7))
#                 title = 'Power-Stats-DA for '+su+' '+conds[0]+' vs '+conds[1]+' '+str(elec)+' '+str(elec_label)+' ('+str(elec_num)+')'
#                 fig.suptitle(title, fontsize=12)

#                 # Plot the POW + STATS
#                 plt.subplot(211)        
#                 BorderPlot(time, x, y=y, kind='sem', alpha=0.2, color=color_codes,linewidth=2, 
#                            ncol=1, xlabel='Time (s)',ylabel = r'Power', legend=conds)
#                 rmaxis(plt.gca(), ['right', 'top'])
#                 addLines(plt.gca(), vLines=[0], vColor=['darkgray'], vWidth=[2])
#                 plt.legend(loc=0, handletextpad=0.1, frameon=False)
#                 plt.gca().yaxis.set_major_locator(MaxNLocator(3,integer=True))

#                 # Plot DA for the POW
#                 plt.subplot(212)
#                 BorderPlot(time, auc, color='b', kind='sd',xlabel='Time (s)', ylim=[0.4,1.], ylabel='Decoding accuracy (%)',linewidth=2, alpha=0.3)
#                 rmaxis(plt.gca(), ['right', 'top'])
#                 addLines(plt.gca(), vLines=[0], vColor=['darkgray'], vWidth=[2])
#                 plt.gca().yaxis.set_major_locator(MaxNLocator(3,integer=True))
#                 plt.plot(time, th_0_05_perm*np.ones(len(time)), '--', color='r', linewidth=2)
#                 #plt.plot(times_plot, th_0_01_perm*np.ones(len(times_plot)), '--', color='orange', linewidth=2)

                #Save plots
                np.save(name_auc, auc)
                np.save(name_perm, perm_scores)
#                 plt.savefig(plot_name, dpi=300, bbox_inches='tight')
#                 plt.clf()
#                 plt.close() 
                del X, auc, pow_data_elec
    del pow_list


# #### ML for the TPSim - ST analysis

# In[16]:


from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedKFold as SKFold
from sklearn.metrics import roc_auc_score
from numpy.random import permutation

conds,phases, subjects = ['low','high'],['odor'],['SEMC','VACJ','FERJ','LEFC','PIRJ','CHAF']#,'MICP','VACJ','SEMC','LEFC','PIRJ']
freq_names = ['2_theta','3_alpha','4_beta','5_gamma1','6_gamma2']
color_codes = ['darkorange','blue']
st = study('Olfacto')
nperm = 1000

for su,phase in product(subjects,phases):
    path_pow = path.join(st.path, 'feature/TPSim_Enc_Ret_By_Odor/TPS_R_p_by_cond/')
    save_path = path.join(st.path, 'classified/Classif_TPSim_Encoding_Retrieval_Odor_1000perm/')

    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(path.join(path_pow, 'TPS_spear_'+su+'_'+phase+'_'+conds[0]+'.npz'))
    names, channels = mat0['label'], mat0['channel']
    #print(mat0['xpow'].shape,time.shape,time[17:42])
    pow_list.append(mat0['TPS']) #nfreq, nelecs, ntrials
    nelecs = mat0['TPS'].shape[1]
    mat1 = np.load(path.join(path_pow, 'TPS_spear_'+su+'_'+phase+'_'+conds[1]+'.npz'))
    pow_list.append(mat1['TPS'])
    print (su, 'power shape: ', [pow.shape for pow in pow_list])
    
    # =========================== Select Power for 1 elec 1 freq =================================                 
    iterator = range(nelecs)
    for elec_num in iterator:
        for freq in range(len(freq_names)):
            elec, elec_label, freq_name = channels[elec_num], names[elec_num], freq_names[freq]
            print ('elec ', elec, 'elec_label ', elec_label)
            #Filenames to save
            name_auc = (save_path+freq_name+'/'+su +'_auc_'+conds[0]+'_'+conds[1]+'_'+str(elec_label)+'_('+str(elec_num)+').npy')
            name_perm = (save_path+freq_name+'/'+su +'_perm_'+str(elec_label)+'_('+str(elec_num)+').npy')
                        
            if not path.exists(name_auc):
                #print(su,phase,elec_num,freq,'already computed')
            #else:
                print('--» processing',su, 'elec', elec_num,'/',nelecs, 'freq',freq)
                pow_data_elec = []
                for i,power in enumerate(pow_list):
                    pow_data_elec.append(power[freq,elec_num][np.newaxis].swapaxes(0,1))
                print('mean TPSim for low & high', [np.mean(power) for power in pow_data_elec])
                nwin = power.shape[1]

        # =============================  Classification Computation ============================================================           
                # create a data matrix, concatenate along the trial dimension
                x = np.concatenate(pow_data_elec, axis=0)
                print ('Size of the concatenated data: ', x.shape, 'Number time windows : ', x.shape[1])
                y = np.hstack([np.array([i]*len(power)) for i, power in enumerate(pow_data_elec)])
                print ('Size of label for classif: ', len(y))

                auc = np.array([])
                for t in range(x.shape[1]):
                    X = x[:,t]
                    X = X.reshape(-1, 1)
                    score_rep = []
                    for i in range(10):
                        k = 5
                        skf = SKFold(n_splits=k, random_state=None, shuffle=True)
                        skf.get_n_splits(X, y)
                        score_cv = []
                        for train_index, test_index in skf.split(X, y):
                            clf = LDA()
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
                        for train_index, test_index in skf.split(X, y_perm):
                            clf = LDA()
                            X_train, X_test = X[train_index], X[test_index]
                            y_train, y_test = y_perm[train_index], y_perm[test_index]
                            clf.fit(X=X_train, y=y_train)
                            y_pred = clf.predict(X_test)
                            score_cv.append(roc_auc_score(y_test,y_pred,average='weighted'))
                        perm_rep.append(np.mean(score_cv))
                    perm_rep = np.asarray(perm_rep).reshape(1,len(perm_rep))
                    perm_scores = np.vstack((perm_scores, perm_rep)) if np.size(perm_scores) else perm_rep
                perm_scores = np.swapaxes(perm_scores,0,1)           
                th_0_05_perm = perm_pvalue2level(perm_scores, p=0.05, maxst=True)
                th_0_01_perm = perm_pvalue2level(perm_scores, p=0.01, maxst=True)
                th_0_001_perm = perm_pvalue2level(perm_scores, p=0.001, maxst=True)
                print('th_perm 05: ', th_0_05_perm[0], '01',th_0_01_perm[0],
                      '001',th_0_001_perm[0], 'auc_mean', np.mean(auc))
                np.save(name_auc, auc)
                np.save(name_perm, perm_scores)
                del X, auc, pow_data_elec
    del pow_list

