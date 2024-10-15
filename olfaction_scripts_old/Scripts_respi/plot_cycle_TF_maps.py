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
from joblib import delayed, Parallel
from Extract_Power_Cond import reshape_to_one_cycle, cycle_no_deform_TF

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
path2save = path_study+'cycleTF_sd/'

subjects = ['LEFC','CHAF','VACJ','SEMC','FERJ','MICP','PIRJ']
sessions = ['E1','E2']
conds = ['odorall','no_odor']
deform = [False,True]

def compute_cycleTF_by_cond(su,sessions,cond,elec_id,deform):
    sess0,sess1 = sessions[0], sessions[1]
    print('deform',cond,deform)
    if deform == True:
        cycle_map0 = all_cycleTF['/{}/{}/{}'.format(su,sess0,elec_id)]
        cycle_map1 = all_cycleTF['/{}/{}/{}'.format(su,sess1,elec_id)]
        print('deform shape', cycle_map0.shape, cycle_map1.shape)
    if deform == False:
        TFmap0 = all_TF['/{}/{}/{}'.format(su,sess0,elec_id)]
        TFmap1 = all_TF['/{}/{}/{}'.format(su,sess1,elec_id)]
        cycle_map0 = cycle_no_deform_TF(TFmap0,su,sess0,nb_point_by_cycle,deform)
        cycle_map1 = cycle_no_deform_TF(TFmap1,su,sess1,nb_point_by_cycle,deform)
        print('no deform shape', cycle_map0.shape, cycle_map1.shape)
  
    #combine all session together
    all_cycle_map0,freqs = reshape_to_one_cycle(cycle_map0,nb_point_by_cycle,deform)
    all_cycle_map1,_ = reshape_to_one_cycle(cycle_map1,nb_point_by_cycle,deform)
    print(all_cycle_map0.shape,all_cycle_map1.shape)
    all_cycle_map = pd.concat((all_cycle_map0,all_cycle_map1), axis=0)
    print(all_cycle_map.shape)
            
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
    print 'size bool', cycle_keep.shape, sum(cycle_keep)

    #select cycles for the condition
    some_cycle_map = all_cycle_map.loc[cycle_keep.astype(bool),:,:]
    n = int(np.round(sum(cycle_keep),0))
    tot = all_cycle_map.shape[0]
    print('n_cycles',n,'/',tot,'for ', su, sess0, cond)
    mean = some_cycle_map.mean(axis=0)
    m_std = some_cycle_map.std(axis=0)
    #m = np.asarray(m_std)
    m = np.asarray(mean)
    df_max = np.max(np.asarray(mean).T)
    clim = (0, df_max)
    return m,n,freqs,clim

def plot_and_savefig_cycle_maps_all():
    for su in subjects:
        sigs = np.load(path_study+'database/'+su+'_E1_sigs_bipo3_phys.npz')
        n_elecs, channels, labels = sigs['x'].shape[0], sigs['channels'], sigs['labels']
        print(n_elecs,su)
        for elec in range(n_elecs):
            elec_id, channel, label = 'elec_'+str(elec), channels[elec], labels[elec]
            plotname = path2save+'{}_{}_{}_np_pts={}_cyclesTF.png'.format(su,'E',channel,nb_point_by_cycle)
            if os.path.isfile(plotname):
                continue
            else:
                map_cond0,n0,freqs,clim0 = compute_cycleTF_by_cond(su,sessions,conds[0],elec_id,deform=True)
                map_cond1,n1,_,_ = compute_cycleTF_by_cond(su,sessions,conds[1],elec_id,deform=True)
                map_cond2,n2,_,clim1 = compute_cycleTF_by_cond(su,sessions,conds[0],elec_id,deform=False)
                map_cond3,n3,_,_ = compute_cycleTF_by_cond(su,sessions,conds[1],elec_id,deform=False)
                map_cond0, map_cond1 = np.asarray(map_cond0), np.asarray(map_cond1)
                map_cond2, map_cond3 = np.asarray(map_cond2), np.asarray(map_cond3)
                #print(map_cond1.T.shape, map_cond1.T[0:20,:].shape)
                
                #Plot all cycle maps on the one figure
                fig, axs = plt.subplots(ncols = 4, figsize = (14,7))
                axs[0].imshow(map_cond0.T, cmap='viridis',  interpolation='bicubic', clim=clim0,
                        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[0].axvline(inspi_ratio, ls = '--', color = 'w')
                axs[0].set_title(conds[0]+' ('+str(n0)+' cycles) - Locked')

                axs[1].imshow(map_cond1.T, cmap='viridis',  interpolation='bicubic', clim=clim0,
                        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[1].axvline(inspi_ratio, ls = '--', color = 'w')
                axs[1].set_title(conds[1]+' ('+str(n1)+' cycles) - Locked')

                axs[2].imshow(map_cond2.T, cmap='viridis',  interpolation='bicubic', clim=clim0,
                       aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[2].set_title(conds[0]+' ('+str(n0)+' cycles) - Unlocked')

                im = axs[3].imshow(map_cond3.T, cmap='viridis',  interpolation='bicubic', clim=clim0,
                      aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[3].set_title(conds[1]+' ('+str(n1)+' cycles) - Unlocked')
                fig.colorbar(im, ax=axs[3], shrink=0.5)
                #cbar = fig.colorbar(plt.gca(), ticks=[-1, 0, 1], orientation='vertical')

                fig.suptitle('{} {} {} {}'.format(su, 'E', channel, label))
                fig.savefig(plotname)
                plt.clf()
                plt.close()

    
if __name__ =='__main__':
   
    plot_and_savefig_cycle_maps_all()



