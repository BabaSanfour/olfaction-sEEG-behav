# -*- coding: utf-8 -*-
import respirationtools
import timefreqtools
import time
import pandas as pd
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import delayed, Parallel

################################ TF Parameters ###########################################################
t_start = 0.
f_start = 0.
f_stop = 120.
nb_freq = 60
f0 = 2.5
norm = 1.
nb_point_by_cycle = 40
inspi_ratio = .4
##########################################################################################################

path_study = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/respi/'
all_TF_name = 'cache_TF/all_TF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_np_pts={}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_TF = pd.HDFStore(path_study+all_TF_name)
all_cycleTF_name = 'cache_TF/all_cycleTF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_np_pts={}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_cycleTF = pd.HDFStore(path_study+all_cycleTF_name)
path_bool = path_study+'cycles_by_cond/'

subjects = ['SEMC','CHAF','VACJ','LEFC','FERJ','MICP','PIRJ']
sessions = ['E1','E2']
conds = ['odorall','no_odor']

def precompute_and_save_cyclefreq(su, sess):
    #Load cycles times
    df_cycles = pd.read_excel(path_study+'cycles_by_cond/All_cycles_features_E.xls',sheetname=su+'_'+sess)
    cycle_times = df_cycles[['insp_time','exp_time']].values
    print(cycle_times.shape)  # (nb_cycle,2)

    #Load signals
    mat = np.load(path_study+'database/{}_{}_sigs_bipo3_phys.npz'.format(su,sess))
    sigs,channels,labels, sf = mat['x'],mat['channels'],mat['Mai_RL'], mat['sf']
    n_elec = sigs.shape[0]
    for elec in range(n_elec):
        print '--> precompute :', su, sess, 'elec',elec
        sig = sigs[elec,:].squeeze()
        
        #Compute TF and save
        wt, times, freqs, sr2 = timefreqtools.compute_timefreq(sig, sf, f_start, f_stop, 
            nb_freq=nb_freq,f0=f0, returns='all', t_start=t_start, min_sampling_rate=f_stop*2.2,
            normalisation=1)
        tfr = pd.DataFrame(np.abs(wt).astype('float32'), index = times, columns=freqs)
        print(tfr.shape)
        
        path = '/{}/{}/{}'.format(su,sess,'elec_'+str(elec))
        all_TF[path] = tfr

        #Compute cycle TF and save
        times = tfr.index.values
        freqs = tfr.columns.values
        
        t1 = time.clock()
        clipped_times, times_to_cycles, cycles, cycle_points, deformed_data = respirationtools.deform_to_cycle_template(tfr.values, 
            times, cycle_times,nb_point_by_cycle=nb_point_by_cycle, inspi_ratio=inspi_ratio)
        t2 = time.clock()
        print('in', t2-t1, 's')
        
        cycle_map = pd.DataFrame(deformed_data.astype('float32'), index=cycle_points, 
            columns=tfr.columns)

        all_cycleTF[path] = cycle_map
        all_TF.flush()
        all_cycleTF.flush()

def precompute_save_all():
    #Parallel(n_jobs=-1)(delayed(precompute_and_save_cyclefreq)(su, sess) for su,sess in product(subjects,sessions))
    for su,sess in product(subjects,sessions):
        precompute_and_save_cyclefreq(su, sess)
        
if __name__ =='__main__':
   
    precompute_save_all()




