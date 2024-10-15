# -*- coding: utf-8 -*-
import respirationtools
import timefreqtools
import time
import pandas as pd
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import os
print respirationtools.__file__

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
all_TF_name = 'cache_TF/tests/all_TF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_np_pts={}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_TF = pd.HDFStore(path_study+all_TF_name)
all_cycleTF_name = 'cache_TF/tests/all_cycleTF_fstart={}_fstop={}_nbfreq={}_f0={}_norm={}_np_pts={}.h5'.format(f_start,f_stop,nb_freq,f0,norm,nb_point_by_cycle)
all_cycleTF = pd.HDFStore(path_study+all_cycleTF_name)
path_bool = path_study+'cycles_by_cond/'
path2save = path_study+'cycles_TF/'

subjects = ['SEMC','CHAF','VACJ','LEFC','FERJ','MICP','PIRJ']
sessions = ['E1','E2']
conds = ['odorall','no_odor']

def precompute_and_save_cyclefreq():
    for su, sess in product(subjects, sessions):
        print(su, sess)

        #Load cycles times
        df_cycles = pd.read_excel(path_study+'cycles_by_cond/All_cycles_features_E.xls',sheetname=su+'_'+sess)
        cycle_times = df_cycles[['insp_time','exp_time']].values
        print(cycle_times.shape)  # (nb_cycle,2)
        

        #Load signals
        mat = np.load(path_study+'database/{}_{}_sigs_bipo3_phys.npz'.format(su,sess))
        sigs,channels,labels, sf = mat['x'],mat['channels'],mat['Mai_RL'], mat['sf']
        n_elec = sigs.shape[0]
    
        for elec in [48]:#range(n_elec):
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
            print('in', t2-t1, 's',cycles,times_to_cycles)
            
            cycle_map = pd.DataFrame(deformed_data.astype('float32'), index=cycle_points, 
                columns=tfr.columns)

            all_cycleTF[path] = cycle_map
            all_cycleTF.flush()


def reshape_to_one_cycle(cycle_map):
    cycle_points =  cycle_map.index.values
    #~ print cycle_points[0], cycle_points[-1]
    
    first_cycle = int(cycle_points[0])
    last_cyle= int(np.ceil(cycle_points[-1]))
    #~ print first_cycle, last_cyle
    cycles = np.arange(first_cycle, last_cyle, dtype='int64')
    
    freqs = cycle_map.columns.values

    onecycle_data = cycle_map.values.reshape(-1, nb_point_by_cycle,  cycle_map.shape[1])
    
    return pd.Panel(onecycle_data, items=cycles, 
                                major_axis=np.arange(nb_point_by_cycle, dtype=float)/nb_point_by_cycle,
                                minor_axis=cycle_map.columns)


def compute_cycleTF_by_cond(su,sessions,cond,elec_id,plot=False):
    sess0,sess1 = sessions[0], sessions[1]
    cycle_map0 = all_cycleTF['/{}/{}/{}'.format(su,sess0,elec_id)]
    cycle_map1 = all_cycleTF['/{}/{}/{}'.format(su,sess1,elec_id)]

    cycle_points0 =  cycle_map0.index.values[:-1]
    cycle_points1 = cycle_map1.index.values[:-1]

    first_cycle0 = int(cycle_points0[0])
    last_cyle0= int(np.round(cycle_points0[-1]))
    first_cycle1 = int(cycle_points1[0])
    last_cyle1= int(np.round(cycle_points1[-1]))
        
    freqs = cycle_map0.columns.values
    
    #combine all session together
    all_cycle_map0 = reshape_to_one_cycle(cycle_map0)
    all_cycle_map1 = reshape_to_one_cycle(cycle_map1)
    print(all_cycle_map0.shape,all_cycle_map1.shape)
    all_cycle_map = pd.concat((all_cycle_map0,all_cycle_map1), axis=0)
    print(all_cycle_map.shape)
            
    if cond == 'all':
        cycle_keep = np.ones(shape=all_cycle_map.shape[0])
    if su == 'VACJ':
        cycle_keep0 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess0,cond))[0:-1]
        cycle_keep1 = np.load(path_bool+ '{}_{}_{}.npy'.format(su,sess1,cond))[16:-1]
        cycle_keep = np.concatenate((cycle_keep0, cycle_keep1),axis=0)
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
    m= some_cycle_map.mean(axis=0)
    m_std = some_cycle_map.std(axis=0)
    m = np.asarray(m)
    df_max0 = np.max(m.T[0:20,:])
    df_max1 = np.max(m.T[30:-1,:])
    clim0, clim1 = (0, df_max0), (0, df_max1)
    
    # m = np.array(m) / np.array(m).mean(1).reshape(-1, 1)

    if plot == True:
        fig, ax = plt.subplots()
        cax = ax.imshow(m.T, cmap='viridis',  interpolation='bicubic', clim=[0, 3.],
                         aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
        title_obj = plt.title('n_cycles :'+str(round(n))+'/'+str(tot)+'\n'+su+' E '+cond)
        cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='vertical')
        ax.axvline(inspi_ratio, ls='--', color='white', lw=2)
        
        plt.show()
        #plt.savefig(path2save+su+'_E_'+cond+'_'+elec+'.png')
        plt.clf()      
        plt.close()
    return m,n,freqs,clim0, clim1

def diff_cycle_maps(su,sessions,conds,elec_id,plot=False):
    cond0, cond1 = conds[0], conds[1]
    map_cond0, n0,freqs,_ = compute_cycleTF_by_cond(su,sessions,cond0,elec_id=elec_id)
    map_cond1, n1,freqs,_ = compute_cycleTF_by_cond(su,sessions,cond1,elec_id=elec_id)

    diff_m = map_cond0-map_cond1#(map_cond0/n0) - (map_cond1/n1)
    max_diff = np.abs(np.array(diff_m)).max()
    clim = (-max_diff, max_diff)
    print(diff_m.shape, clim)

    
    if plot == True:
        fig, ax = plt.subplots()
        cax = ax.imshow(diff_m.T, cmap='viridis',  interpolation='bicubic', clim=[0,0.05],
                         aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
        title_obj = plt.title(cond0[:4]+'('+str(n0)+') - '+cond1+'('+str(n1)+') // '+su+' in '+elec_id)
        #cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='vertical')
        ax.axvline(inspi_ratio, ls = '--', color = 'w')
        plt.savefig(path2save+su+'_E_'+cond0[:4]+'_vs_'+cond1+'_'+elec_id+'.png')
        plt.clf()
        plt.close()

    return diff_m, clim

def plot_and_savefig_diff_cycle_maps_all():
    for su in subjects:
        sigs = np.load(path_study+'database/'+su+'_E1_sigs_bipo3_phys.npz')
        n_elecs, channels, labels = sigs['x'].shape[0], sigs['channels'], sigs['labels']
        print(n_elecs,su)
        for elec in range(n_elecs):
            elec_id, channel, label = 'elec_'+str(elec), channels[elec], labels[elec]
            plotname = path2save+'{}_{}_{}_{}_cyclesTF_low.png'.format(su,'E',channel,label)
            if os.path.isfile(plotname):
                continue
            else:
                map_cond0,n0,freqs,clim0,clim1 = compute_cycleTF_by_cond(su,sessions,conds[0],elec_id)
                map_cond1,n1,_,_,_ = compute_cycleTF_by_cond(su,sessions,conds[1],elec_id)
                #map_cond2,n2,_,_,_ = compute_cycleTF_by_cond(su,sessions,conds[2],elec_id)
                map_cond0, map_cond1 = np.asarray(map_cond0), np.asarray(map_cond1)
                #map_cond2 = np.asarray(map_cond2)
                #print(map_cond1.T.shape, map_cond1.T[0:20,:].shape)
                
                #map_diff0,clim = diff_cycle_maps(su,sessions,conds[:2],elec_id)
                #map_diff1,clim = diff_cycle_maps(su,sessions,[conds[2],conds[1]],elec_id)
                #map_diff2,clim = diff_cycle_maps(su,sessions,[conds[0],conds[2]],elec_id)
                
                #Plot all cycle maps on the one figure
                fig, axs = plt.subplots(ncols = 2, figsize = (8,7))
                axs[0].imshow(map_cond0.T, cmap='viridis',  interpolation='bicubic', clim=clim1,
                        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[0].axvline(inspi_ratio, ls = '--', color = 'w')
                axs[0].set_title(conds[0]+' ('+str(n0)+' cycles) - LOW')

                axs[1].imshow(map_cond1.T, cmap='viridis',  interpolation='bicubic', clim=clim1,
                        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                axs[1].axvline(inspi_ratio, ls = '--', color = 'w')
                axs[1].set_title(conds[1]+' ('+str(n1)+' cycles) - LOW')

                #axs[2].imshow(map_cond0.T[30:-1,:], cmap='viridis',  interpolation='bicubic', clim=clim1,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[30], freqs[-1]))
                #axs[2].axvline(inspi_ratio, ls = '--', color = 'w')
                #axs[2].set_title(conds[0]+' ('+str(n0)+' cycles) - HIGH')

                #axs[3].imshow(map_cond1.T[30:-1,:], cmap='viridis',  interpolation='bicubic', clim=clim1,
                #       aspect = 'auto',origin='lower', extent=(0, 1, freqs[30], freqs[-1]))
                #axs[3].axvline(inspi_ratio, ls = '--', color = 'w')
                #axs[3].set_title(conds[1]+' ('+str(n1)+' cycles) - HIGH')
                
                # axs[4].imshow(map_cond2.T[30:-1,:], cmap='viridis',  interpolation='bicubic', clim=clim1,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[30], freqs[-1]))
                # axs[4].axvline(inspi_ratio, ls = '--', color = 'w')
                # axs[4].set_title(conds[2]+' - '+conds[1])

                # axs[5].imshow(map_cond1.T[30:-1,:], cmap='viridis',  interpolation='bicubic', clim=clim1,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[30], freqs[-1]))
                # axs[5].axvline(inspi_ratio, ls = '--', color = 'w')
                # axs[5].set_title(conds[0]+' - '+conds[2])

                #axs[3].imshow(map_diff0.T, cmap='bwr',  interpolation='bicubic', clim=clim,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                #axs[3].axvline(inspi_ratio, ls = '--', color = 'w')
                #axs[3].set_title(conds[0]+' - '+conds[1])
                
                #axs[4].imshow(map_diff1.T, cmap='bwr',  interpolation='bicubic', clim=clim,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                #axs[4].axvline(inspi_ratio, ls = '--', color = 'w')
                #axs[4].set_title(conds[2]+' - '+conds[1])

                #axs[5].imshow(map_diff2.T, cmap='bwr',  interpolation='bicubic', clim=clim,
                #        aspect = 'auto',origin='lower', extent=(0, 1, freqs[0], freqs[-1]))
                #axs[5].axvline(inspi_ratio, ls = '--', color = 'w')
                #axs[5].set_title(conds[0]+' - '+conds[2])
             
                fig.suptitle('{} {} {} {}'.format(su, 'E', channel, label))
                #plt.show()
                fig.savefig(plotname)
                plt.clf()
                plt.close()

if __name__ =='__main__':
   
    precompute_and_save_cyclefreq()
    #test_plot_by_cond()
    #compute_cycleTF_by_cond('LEFC',['E1','E2'],'all','elec_0',plot=True)
    #diff_cycle_maps(su='LEFC',sessions=['E1','E2'],conds=['odor0','no_odor'],elec_id='elec_0',plot=False)
    #plot_and_savefig_diff_cycle_maps_all()



