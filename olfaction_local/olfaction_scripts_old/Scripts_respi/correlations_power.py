# -*- coding: utf-8 -*-
import respirationtools
import timefreqtools
import time
import pandas as pd
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import os
#print respirationtools.__file__

################################ TF Parameters ###########################################################
f_start, f_stop = 0., 120.
norm, nb_freq, f0 = 1., 60, 2.5
nb_point_by_cycle = 120
inspi_ratio = .4
##########################################################################################################
path_study = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/respi/'
all_cycleTF_name = 'cache_TF/all_cycleTF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_{}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_cycleTF = pd.HDFStore(path_study+all_cycleTF_name)
all_TF_name = 'cache_TF/all_TF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_{}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_TF = pd.HDFStore(path_study+all_TF_name)
path_bool = path_study+'cycles_by_cond/'
path2save = path_study+'Power_reshape/test/'
##########################################################################################################

subjects = ['LEFC']#,'CHAF','VACJ','SEMC','FERJ','MICP','PIRJ']
sessions = ['E1','E2']
bands = [(0,1),(1,2),(2,4),(4,8),(8,18),(33,59)]
freqnames = ['VLFC', 'delta','theta','alpha','beta','gamma2']
conds = ['odorall','no_odor']
deform = True

def reshape_to_one_cycle(cycle_map):
    cycle_points =  cycle_map.index.values
    
    first_cycle = int(cycle_points[0])
    last_cyle= int(np.ceil(cycle_points[-1]))
    #print(first_cycle,last_cyle)
    cycles = np.arange(first_cycle, last_cyle, dtype='int64')
    
    freqs = cycle_map.columns.values
    onecycle_data = cycle_map.values.reshape(-1, nb_point_by_cycle,  cycle_map.shape[1])
    print 'reshape', onecycle_data.shape
    
    return pd.Panel(onecycle_data, items=cycles, 
                        major_axis=np.arange(nb_point_by_cycle, dtype=float)/nb_point_by_cycle,
                        minor_axis=cycle_map.columns)

def cycle_TF(TF_map,su,sess):
    #Compute fake cycle TF
    times = TF_map.index.values
    freqs = TF_map.columns
    
    df_cycles = pd.read_excel(path_study+'cycles_by_cond/All_cycles_features_E.xls',sheetname=su+'_'+sess)
    start_times = df_cycles[['insp_time']].values
    end_times = start_times + 1.3
    cycle_times = np.concatenate((start_times,end_times), axis=1)
    #print(cycle_times.shape, cycle_times[:3])

    clipped_times, times_to_cycles, cycles, cycle_points, deformed_data = respirationtools.deform_to_cycle_template(TF_map.values,
    times, cycle_times,nb_point_by_cycle=nb_point_by_cycle, inspi_ratio=inspi_ratio)
    cycle_map = pd.DataFrame(deformed_data.astype('float32'), index=cycle_points, 
                columns=freqs)
    #print cycle_map.shape
    return cycle_map

def export_power_by_cond(su,sessions,cond,elec_id,deform=True):
    sess0,sess1 = sessions[0], sessions[1]
    
    if deform == True:
        cycle_map0 = all_cycleTF['/{}/{}/{}'.format(su,sess0,elec_id)]
        cycle_map1 = all_cycleTF['/{}/{}/{}'.format(su,sess1,elec_id)]
    if deform == False:
        TFmap0 = all_TF['/{}/{}/{}'.format(su,sess0,elec_id)]
        TFmap1 = all_TF['/{}/{}/{}'.format(su,sess1,elec_id)]
        cycle_map0 = cycle_TF(TFmap0,su,sess0)
        cycle_map1 = cycle_TF(TFmap1,su,sess1)

    #combine all sessions together
    all_cycle_map0 = reshape_to_one_cycle(cycle_map0)
    all_cycle_map1 = reshape_to_one_cycle(cycle_map1)
    all_cycle_map = pd.concat((all_cycle_map0,all_cycle_map1), axis=0)
    print(all_cycle_map.shape)

    #select cycles for a specific cond
    if cond == 'all':
        cycle_keep = np.ones(shape=all_cycle_map.shape[0])
    if su == 'VACJ':
        cycle_keep0 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess0,cond))[0:-1]
        cycle_keep1 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess1,cond))[16:-1]
        cycle_keep = np.concatenate((cycle_keep0, cycle_keep1),axis=0)
        if cycle_keep.shape != all_cycle_map.shape[0]:
            cycle_keep = cycle_keep[:-1]
    else:
        cycle_keep0 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess0,cond))[0:-1]
        cycle_keep1 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess1,cond))[0:-1]
        cycle_keep = np.concatenate((cycle_keep0, cycle_keep1),axis=0)
        if cycle_keep.shape != all_cycle_map.shape[0]:
            cycle_keep = cycle_keep[:-1]
    print 'size bool', cycle_keep.shape, sum(cycle_keep), 'for',cond
    some_cycle_map = all_cycle_map.loc[cycle_keep.astype(bool),:,:]
    some_cycle_map = np.asarray(some_cycle_map)
        
    #extract power values by frequency bands
    elec_pow = np.array([])
    for b in range(len(bands)):
        f_pow = np.mean(some_cycle_map[:,:,bands[b]],axis=2)[np.newaxis]
        elec_pow = np.vstack((elec_pow,f_pow)) if np.size(elec_pow) else f_pow
    return elec_pow

def export_all_power():
    for su, cond in product(subjects,conds):
        if deform == True:
            save_x = path2save+ su+'_E_'+cond+'_pow_reshape_120.npz'
        else:
            save_x = path2save+ su+'_E_'+cond+'_pow_120.npz'
        if os.path.isfile(save_x):
            print(su, cond, 'already computed')
        else:
            #Create a dict for all parameters
            kwargs = {}
            kwargs['f'] = [[0.,2],[2,4],[4,8],[8,16],[16,36],[66,120]]
            kwargs['fname'] = ['VLFC','delta','theta','alpha','beta','gamma2']
            sigs = np.load(path_study+'database/'+su+'_E1_sigs_bipo3_phys.npz')
            n_elecs, kwargs['channels'], kwargs['labels'] = sigs['x'].shape[0], sigs['channels'], sigs['Mai_RL']
            kwargs['xyz'], kwargs['aal'], kwargs['BA'] = sigs['xyz'], sigs['aal'], sigs['BA']
            print('--> processing', su, cond, 'n_elecs', n_elecs)

            #Compute power features for all freqs, all elecs
            su_pow = np.array([])
            for elec in [72]:#range(n_elecs):
                elec_id = 'elec_'+str(elec)
                power = export_power_by_cond(su,sessions,cond,elec_id,deform=deform)         
                su_pow = np.vstack((su_pow,power[np.newaxis])) if np.size(su_pow) else power[np.newaxis]
            kwargs['pow'] = su_pow
            np.savez(save_x, **kwargs)

if __name__ =='__main__':
   
    export_all_power()



