
# coding: utf-8

# ## Import Libraries

# In[ ]:

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

from os import path
from mne.stats import *
from mne.baseline import rescale
from mne.filter import filter_data
import time


# ## User variables

# In[2]:

# where to find data
st = study('Olfacto')
score = 'Epi' #'Rec'
if score == 'Epi':
    path_elec = path.join(st.path, 'database/TS_E_all_by_odor_th40_art400_30_250_5s_Good_Bad_EpiScore/')
    path_pow = path.join(st.path, 'feature/7_Power_E1E2_Odor_Good_Bad_700_100_EpiScore/')
    save_path = path.join(st.path, 'classified/1_Classif_Power_EpiScore_all_electrodes_win700_step100/')

# ANALYSIS PARAMETERS
classif = 'lda'
nfreq = 5
alpha = 0.05
minsucc = 3 #nb of continuous samples to be significant


# ## Power Decoding - Good Bad Odors Encoding

# In[3]:

test = False

if test == True:
    n_elec = {'PIRJ' :1}
    subjects = ['PIRJ']
else :
    subjects = ['MICP','VACJ','SEMC','PIRJ','LEFC','CHAF'] 
    n_elec = {
    'CHAF' : 81,
    'VACJ' : 91, 
    'SEMC' : 81,
    'PIRJ' : 62,
    'LEFC' : 160,
    'MICP' : 79,
        }

for su in subjects:
    #load power files (nfreq, nelec, nwin, ntrial)
    bad_data = np.load(path.join(path_pow, su+'_concat_odor_bad_bipo_new_power.npz'))['xpow']
    good_data = np.load(path.join(path_pow, su+'_concat_odor_good_bipo_new_power.npz'))['xpow']
    print (su, 'bad shape: ', bad_data.shape, 'good shape: ', good_data.shape)
    
# ==========================  BALANCED CONDITIONS - Bootstrap  =====================================
    if bad_data.shape[3] > good_data.shape[3]:
        bad_stat = bad_data[:,:,:,np.random.randint(bad_data.shape[3], size=good_data.shape[3])]
        good_stat = good_data
    elif bad_data.shape[3] < good_data.shape[3]:
        bad_stat = bad_data
        good_stat = good_data[:,:,:,np.random.randint(good_data.shape[3], size=bad_data.shape[3])]
    else:
        bad_stat, good_stat = bad_data, good_data
    ntrials = bad_stat.shape[3]
    print ('balanced data : ', bad_stat.shape, good_stat.shape)
                 
# =========================== SELECT Power for 1 elec 1 freq =================================                 
    for elec_num in range(n_elec[su]):
        for freq in range(nfreq):
            # load power files for 1 elec // 1 freq // Bad-Good conditions
            bad_data_elec = bad_stat[freq,elec_num].swapaxes(0,1)
            good_data_elec = good_stat[freq,elec_num].swapaxes(0,1)
            print ('data elec ', bad_data_elec.shape, good_data_elec.shape)
            nwin = good_data.shape[1]
            elec = np.load(path.join(path_elec, su+'_concat_odor_bad_bipo_new.npz'))['channel'][elec_num]
            elec_label = np.load(path.join(path_elec, su+'_concat_odor_bad_bipo_new.npz'))['label'][elec_num]
            freq_name = np.load(path.join(path_pow, su+'_concat_odor_bad_bipo_new_power.npz'))['fname'][freq]
            print ('elec ', elec, 'elec_label ', elec_label)

# ===========================  STATISTICS  =====================================
            # Permutations and t test of the data
            bad_perm, good_perm = perm_swap(bad_data_elec, good_data_elec, n_perm=1000, axis=0)
            bad_perm, good_perm = np.swapaxes(bad_perm,0,1), np.swapaxes(good_perm,0,1)
            print('data permuted', bad_perm.shape, good_perm.shape)
            Tperm, _ = ttest_ind(bad_perm, good_perm, equal_var=False)
            print('T perm', Tperm.shape)
            thr_0_5_stat = [-perm_pvalue2level(Tperm, p=0.05, maxst=True)[0],perm_pvalue2level(Tperm, p=0.05, maxst=True)[0]]
            thr_0_1_stat = [-perm_pvalue2level(Tperm, p=0.01, maxst=True)[0],perm_pvalue2level(Tperm, p=0.01, maxst=True)[0]]
            thr_0_0_1_stat = [-perm_pvalue2level(Tperm, p=0.001, maxst=True)[0],perm_pvalue2level(Tperm, p=0.001, maxst=True)[0]]
            print('treshold stats', thr_0_5_stat,thr_0_1_stat,thr_0_0_1_stat)
            T0, _  = ttest_ind(bad_data_elec, good_data_elec, equal_var=False)
            print('Obs stats',T0.shape, T0.max(), T0.min())

            # Create the pvalue vector to plot
            pvals = []
            for i in range(T0.shape[0]):
                if T0[i] < thr_0_0_1_stat[0] or T0[i] > thr_0_0_1_stat[1]:
                    pval = pvals.append(0.0009)
                elif T0[i] < thr_0_1_stat[0] or T0[i] > thr_0_1_stat[1]:
                    pval = pvals.append(0.009)
                elif T0[i] < thr_0_5_stat[0] or T0[i] > thr_0_5_stat[1]:
                    pval = pvals.append(0.04)
                else:
                    pval = pvals.append(1)
            print (pvals)
                
# =============================  CLASSIFICATION COMPUTATION ============================================================           
            #create a data matrix, concatenate along the trial dimension
            bad_good = np.concatenate((bad_data_elec, good_data_elec), axis=0)
            print ('Size of the concatenated data: ', bad_good.shape, 'Number time windows : ', bad_good.shape[1])
            #create label vector (0 for rest and 1 for odor)
            y = [0]*bad_data_elec.shape[0] + [1]*good_data_elec.shape[0]
            print ('Size of label for classif: ', len(y))
            # Define a cross validation:
            cv = defCv(y, n_folds=10, cvtype='skfold', rep=10)
            # Define classifier technique
            clf = defClf(y=y, clf=classif)#,n_tree=200, random_state=100)
            #Classify rest and odor
            cl = classify(y, clf=clf, cvtype=cv)
            # Evaluate the classifier on data:
            da,pvalues,daperm = cl.fit(bad_good, n_perm=1000,method='full_rnd', mf=False)
            #print(pvalues.shape, pvalues.min(), pvalues.max())
            th_0_05_perm = perm_pvalue2level(daperm, p=0.05, maxst=True)
            th_0_01_perm = perm_pvalue2level(daperm, p=0.01, maxst=True)
            th_0_001_perm = perm_pvalue2level(daperm, p=0.001, maxst=True)
            print('th_perm : ', th_0_05_perm[0], th_0_01_perm[0], th_0_001_perm[0])

    # ============================== PLOT POWER ANALYSIS + STATS & DECODING ACCURACY ===================================================
            # plot and figure parameters
            xfmt = ScalarFormatter(useMathText=True)
            xfmt.set_powerlimits((0,3))
            fig = plt.figure(1,figsize=(7,7))
            title = 'Power-Stats-DA for '+su+' Good/Bad '+str(elec)+' '+str(elec_label)+' ('+str(elec_num)+') ntrials:'+str(ntrials)
            fig.suptitle(title, fontsize=12)
            # Time vector to plot power
            step = 3500/ bad_data_elec.shape[1]
            times_plot = np.arange(-500, 3000, step)

            # Plot the POW + STATS
            plt.subplot(211)
            bad_good_to_plot = bad_good * 100
            BorderPlot(times_plot, bad_good_to_plot, y=y, kind='sem', alpha=0.2, color=['b','m'], 
                       linewidth=2, ncol=1, xlabel='Time (ms)',ylabel = r'Power change (%)', legend=['bad','good'])
            addLines(plt.gca(), vLines=[0], vColor=['r'], vWidth=[2], hLines=[0], 
                     hColor=['#000000'], hWidth=[2])
            addPval(plt.gca(), pvals, p=0.05, x=times_plot, y=5, color='orange', lw=2, minsucc=minsucc)
            addPval(plt.gca(), pvals, p=0.01, x=times_plot, y=5, color='r', lw=2,minsucc=minsucc)
            addPval(plt.gca(), pvals, p=0.001, x=times_plot, y=5, color='g', lw=2,minsucc=minsucc)
            rmaxis(plt.gca(), ['right', 'top'])
            plt.legend(loc=0, handletextpad=0.1, frameon=False)
            plt.gca().yaxis.set_major_locator(MaxNLocator(3,integer=True))

            # Plot DA for the POW
            plt.subplot(212)
            BorderPlot(times_plot, da, color='b', kind='sem',xlabel='Time (ms)', 
                       ylim=[da.min()-10,da.max()+10], ylabel='Decoding accuracy (%)',
                       linewidth=2, alpha=0.3)
            rmaxis(plt.gca(), ['right', 'top'])
            addLines(plt.gca(), vLines=[0], vWidth=[2], vColor=['r'], hLines=[50], 
                     hColor=['#000000'], hWidth=[2])
            plt.legend(loc=0, handletextpad=0.1, frameon=False)   
            plt.gca().yaxis.set_major_locator(MaxNLocator(3,integer=True))
            plt.plot(times_plot, th_0_05_perm*np.ones(len(times_plot)), '--', color='orange', linewidth=2)
            plt.plot(times_plot, th_0_01_perm*np.ones(len(times_plot)), '--', color='r', linewidth=2)
            plt.plot(times_plot, th_0_001_perm*np.ones(len(times_plot)), '--', color='g', linewidth=2)
            # Criteria to be significant
            pvals = np.ravel(pvals)
            underp = np.where(pvals < alpha)[0]
            pvsplit = np.split(underp, np.where(np.diff(underp) != 1)[0]+1)
            signif = [True for k in pvsplit if len(k) >= minsucc]
        
            #Save plots and stats
            if len(signif) >=1:
                name_t0 = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/stat/'+su +'_t0_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_tperm = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/stat/'+su +'_tperm_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_pval = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/stat/'+su +'_pvals_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_da = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/da/'+su +'_da_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_daperm = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/da/'+su +'_daperm_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                plot_name = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Significant/fig/'+su +'_Power_'  + score +'_'+str(elec_label)+'_('+str(elec_num)+').png')
            else:
                name_t0 = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/stat/'+su +'_t0_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_tperm = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/stat/'+su +'_tperm_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_pval = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/stat/'+su +'_pvals_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_da = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/da/'+su +'_da_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                name_daperm = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/da/'+su +'_daperm_' + score +'_'+str(elec_label)+'_('+str(elec_num)+').npy')
                plot_name = (save_path+'All_balanced_1_1000perm_DA_stats_sametrials/'+str(freq)+'_'+freq_name+'/Not_Significant/fig/'+su +'_Power_'  + score +'_'+str(elec_label)+'_('+str(elec_num)+').png')

            np.save(name_t0, T0)
            np.save(name_tperm, Tperm)
            np.save(name_pval, pvals)
            np.save(name_da, da)
            np.save(name_daperm, daperm)
            plt.savefig(plot_name, dpi=300, bbox_inches='tight')
            plt.clf()
            plt.close() 
            del bad_data_elec, good_data_elec, bad_perm, good_perm, bad_good, da, pvalues, daperm,
    del bad_data, good_data, bad_stat, good_stat
        

