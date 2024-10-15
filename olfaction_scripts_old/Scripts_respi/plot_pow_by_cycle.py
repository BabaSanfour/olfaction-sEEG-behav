# -*- coding: utf-8 -*-

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
nb_point_by_cycle = 60
inspi_ratio = .4
##########################################################################################################
path_study = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/respi/'
path_pow = path_study+'Power_reshape/'
path2save = path_study+'cycleTF_sd/'
##########################################################################################################

subjects = ['LEFC']
conds = ['odorall','no_odor']
deforms = ['_reshape','']

def plot_pow_by_cycle():
    pow_plot, cond_plot = [],[]
    for su, cond, deform in product(subjects,conds, deforms):
        mat = np.load(path_pow+ '{}_E_{}_pow{}_{}.npz'.format(su,cond,deform,nb_point_by_cycle))
        channels, labels = mat['channels'], mat['labels']
        elec = 0
        elec_id, channel, label = 'elec_'+str(elec), channels[elec], labels[elec]

        #power dimensions // nelecs, nfreqs, ncycles, ntimepoints
        pow_plot.append(mat['pow'][elec,4,:,:])
        cond_plot.append(cond+deform)
    print([x.shape for x in pow_plot],[c for c in cond_plot])
    
    #Plot all power cycles by cond
    plotname = path2save+'{}_{}_{}_np_pts={}_Plot.png'.format(su,'E',channel,nb_point_by_cycle)
    fig, axs = plt.subplots(ncols = 4, figsize = (20,7))
    axs[0].plot(pow_plot[0].T)
    axs[0].axvline(inspi_ratio*nb_point_by_cycle, ls = '--', color = 'black')
    axs[0].set_title('('+str(pow_plot[0].shape)+' cycles) - '+cond_plot[0])

    axs[1].plot(pow_plot[1].T)
    axs[1].set_title('('+str(pow_plot[1].shape)+' cycles) - '+cond_plot[1])

    axs[2].plot(pow_plot[2].T)
    axs[2].axvline(inspi_ratio*nb_point_by_cycle, ls = '--', color = 'black')
    axs[2].set_title('('+str(pow_plot[2].shape)+' cycles) - '+cond_plot[2])

    axs[3].plot(pow_plot[3].T)
    axs[3].set_title('('+str(pow_plot[3].shape)+' cycles) - '+cond_plot[3])

    fig.suptitle('{} {} {} {}'.format(su, 'E', channel, label))
    fig.savefig(plotname)
    plt.clf()
    plt.close()

def plot_single_trial_pow():
    pow_plot, cond_plot, n_cycles = [],[],[]
    for su, cond, deform in product(subjects,conds, deforms):
        mat = np.load(path_pow+ '{}_E_{}_pow{}_{}.npz'.format(su,cond,deform,nb_point_by_cycle))
        channels, labels = mat['channels'], mat['labels']
        elec = 0
        elec_id, channel, label = 'elec_'+str(elec), channels[elec], labels[elec]

        #power dimensions // nelecs, nfreqs, ncycles, ntimepoints
        pow_cond = mat['pow'][elec,4,:,:]
        pow_plot.append(pow_cond)
        cond_plot.append(cond+deform)
        n_cycles.append(pow_cond.shape[0])
    print([x.shape for x in pow_plot],[c for c in cond_plot])
    
    #Plot all power cycles by cond
    plotname = path2save+'{}_{}_{}_np_pts={}_cyclesPlot.png'.format(su,'E',channel,nb_point_by_cycle)
    fig, axs = plt.subplots(ncols = 4, figsize = (20,7))
    axs[0].imshow(pow_plot[0],cmap='jet',aspect = 'auto',origin='lower', extent=(0,120,0,n_cycles[0]),
        clim=(0.1,0.3))
    axs[0].axvline(inspi_ratio*nb_point_by_cycle, ls = '--', color = 'w')
    axs[0].set_title('('+str(pow_plot[0].shape)+' cycles) - '+cond_plot[0])

    axs[1].imshow(pow_plot[1],cmap='jet',aspect = 'auto',origin='lower', extent=(0,120,0,n_cycles[0]),
        clim=(0.1,0.3))
    axs[1].set_title('('+str(pow_plot[1].shape)+' cycles) - '+cond_plot[1])

    axs[2].imshow(pow_plot[2],cmap='jet',aspect = 'auto',origin='lower', extent=(0,120,0,n_cycles[0]),
        clim=(0.1,0.3))
    axs[2].axvline(inspi_ratio*nb_point_by_cycle, ls = '--', color = 'w')
    axs[2].set_title('('+str(pow_plot[2].shape)+' cycles) - '+cond_plot[2])

    axs[3].imshow(pow_plot[3],cmap='jet',aspect = 'auto',origin='lower', extent=(0,120,0,n_cycles[0]),
        clim=(0.1,0.3))
    axs[3].set_title('('+str(pow_plot[3].shape)+' cycles) - '+cond_plot[3])

    fig.suptitle('{} {} {} {}'.format(su, 'E', channel, label))
    fig.savefig(plotname)
    plt.clf()
    plt.close()



if __name__ =='__main__':
   
    plot_pow_by_cycle()
    #plot_single_trial_pow()



