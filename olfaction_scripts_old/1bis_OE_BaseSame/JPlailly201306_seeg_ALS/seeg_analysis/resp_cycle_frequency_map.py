# -*- coding: utf-8 -*-
"""
Etude des cartes temps frequences avec deformation sur cycles respi.
"""
#~ import matplotlib
#~ matplotlib.use('Agg')

import sys, os
sys.path.append('../behavior') 
from connection import *

import quantities as pq
import numpy as np
import scipy.signal
from OpenElectrophy.timefrequency import TimeFreq

import pandas as pd

from matplotlib import pyplot
import time
import neo

import scipy.interpolate
import multiprocessing as mp


from artifact_detection import is_artifact_in_window


import params
from joblib import Memory
memory  = Memory(cachedir=params.joblib_cachedir)


f_start = 0.
f_stop = 120.
deltafreq = .3
f0 = 2.5
sampling_rate_tfr = f_stop*3

clim = [0, 30]

f_start_env = 1.
f_stop_env = 40.

nb_point_by_cycle = 30
inspi_ratio = .4
nb_point_inspi = int(nb_point_by_cycle*inspi_ratio)
nb_point_expi = nb_point_by_cycle - nb_point_inspi
one_cycle = np.linspace(0,1,nb_point_by_cycle)
two_cycle = np.linspace(0,2,nb_point_by_cycle*2)

# mean_cycles_maps_save_dirname = './mean_cycles_maps'
# if not os.path.exists(mean_cycles_maps_save_dirname):
#     os.mkdir(mean_cycles_maps_save_dirname)
    
diff_cycle_maps_save_dirname = './diff_cycles_maps_0-40'
if not os.path.exists(diff_cycle_maps_save_dirname):
    os.mkdir(diff_cycle_maps_save_dirname)

# bytrial_envelop_cycles_save_dirname = './bytrial_envelop_cycles'
# if not os.path.exists(bytrial_envelop_cycles_save_dirname):
#     os.mkdir(bytrial_envelop_cycles_save_dirname)


def plot_map(map, ax, times, clim = None, colorbar = True, cax = None):
    im = ax.imshow(map.transpose(),
                                interpolation='nearest', 
                                extent=(times[0], times[-1], f_start-deltafreq/2., f_stop-deltafreq/2.),
                                origin ='lower' ,
                                aspect = 'auto')
    if clim is not None:
        im.set_clim(clim)
    if colorbar:
        if cax is None:
            ax.figure.colorbar(im)
        else:
            ax.figure.colorbar(im,ax = ax, cax = cax ,orientation=orientation)

def estimate_inspi_ratio():
    all_ratio = []
    for subject in session.query(Subject):
        for run in subject.runs:
            seg = run.segments[0]
            respsig = seg.respirationsignals[0]
            cycle_times = respsig.cycle_times.magnitude
            ratio = (cycle_times[:-1,1]-cycle_times[:-1,0])/(cycle_times[1:,0]-cycle_times[:-1,0])
            all_ratio.append(ratio)
    all_ratio = np.concatenate(all_ratio, axis=0)
    fig, ax = pyplot.subplots()
    y,x  = np.histogram(all_ratio, bins = np.arange(0,1,0.01))
    ax.plot(x[:-1], y)
    ax.axvline(np.mean(all_ratio), ls = '-')
    ax.axvline(np.median(all_ratio), ls = '-')
    pyplot.show()


#~ @memory.cache
def compute_full_cycle_map(anasig_id, respsig_id, f_start, f_stop, deltafreq, f0, nb_point_by_cycle):
    dbinfo = open_db(url, myglobals = locals(), use_global_session = False, 
                            object_number_in_cache = 0,  numpy_storage_engine = 'sqltable', compress = compress,
                            relationship_lazy = 'dynamic', predefined_classes = classes)
    session = dbinfo.Session()
    
    anasig =  session.query(AnalogSignal).get(anasig_id)
    respsig =  session.query(RespirationSignal).get(respsig_id)
    
    cycle_times = respsig.cycle_times.magnitude
    #~ print cycle_times.shape
    nb_cycles = cycle_times.shape[0]-1
    
    neosig = anasig.to_neo()
    #~ neosig.oeinstance = None # this release the db object to avoid conenciton to keep open
    tfr = TimeFreq(neosig,
                                t_start = neosig.t_start,
                                t_stop =neosig.times[-1],
                                f_start = f_start * pq.Hz,
                                f_stop = f_stop* pq.Hz,
                                deltafreq =deltafreq* pq.Hz,
                                f0 = f0,
                                sampling_rate =  sampling_rate_tfr * pq.Hz,
                                optimize_fft = True,
                                use_joblib = False,
                                )
    map = np.abs(tfr.map).astype('float32')
    full_times =  tfr.times.rescale('s').magnitude

    keep = (full_times>=cycle_times[0,0]) & (full_times<cycle_times[-1,0])
    map2 = map[keep, :]
    times = full_times[keep]

    # construct cycle_step
    first_cycle = None
    time_to_cycles = np.zeros(times.shape)*np.nan
    for c in range(nb_cycles):
        #~ #One segment
        #~ mask = (times>=cycle_times[c, 0]) & (times<cycle_times[c+1, 0])
        #~ time_to_cycles[mask] = np.linspace(c, c+1, num = np.sum(mask))
        
        #2 segments : inspi + expi
        if cycle_times[c, 0]<times[0]:
            continue
        if first_cycle is None:
            first_cycle = c
        mask = (times>=cycle_times[c, 0]) & (times<cycle_times[c, 1])
        #~ time_to_cycles[mask] = np.linspace(c, c+inspi_ratio, num = np.sum(mask))
        time_to_cycles[mask] = np.arange(np.sum(mask), dtype = float)/np.sum(mask)*inspi_ratio+c
        mask = (times>=cycle_times[c, 1]) & (times<cycle_times[c+1, 0])
        #~ time_to_cycles[mask] = np.linspace(c+inspi_ratio, c+1, num = np.sum(mask))
        time_to_cycles[mask] = np.arange(np.sum(mask), dtype = float)/np.sum(mask)*(1-inspi_ratio)+c+inspi_ratio
    
    keep = ~np.isnan(time_to_cycles)
    time_to_cycles = time_to_cycles[keep]
    times = times[keep]
    map2 = map2[keep]
    
    interp = scipy.interpolate.interp1d(time_to_cycles, map2, kind='linear', axis=0)
    all_cycles = np.arange(first_cycle, nb_cycles, 1./nb_point_by_cycle)

    if all_cycles[-1]>time_to_cycles[-1]:
        nb_cycles = nb_cycles-1
        all_cycles = np.arange(first_cycle, nb_cycles, 1./nb_point_by_cycle)
    
    cycle_map = interp(all_cycles)
        
    session.close()
    dbinfo.Session.get_bind().dispose()
    return tfr, first_cycle, nb_cycles,  cycle_map, time_to_cycles, times, cycle_times, all_cycles


def sum_cycles(cycle_map, cycles_list, first_cycle):
    shape = list(cycle_map.shape)
    shape[0] = nb_point_by_cycle
    sum_cycle_map = np.zeros(shape)
    for c in cycles_list:
        c2 = c - first_cycle
        sum_cycle_map += cycle_map[c2*nb_point_by_cycle:(c2+1)*nb_point_by_cycle]
    return sum_cycle_map

def mean_cycles(cycle_map, cycles_list, first_cycle):
    return sum_cycles(cycle_map, cycles_list, first_cycle)/len(cycles_list)



def get_cycle_nums(run, anasig, first_cycle, nb_cycles, cycle_times):
    trigs = [trial.triggered_odor_time for trial in run.trials.order_by(text('`index`'))]
    cycle_before=[]
    cycle_after=[]
    for trig in trigs:
	if trig>cycle_times[-1,0]: continue #LEFC R2 en dehors du TRC mais essai OK du coup inclus
        c0 = np.where(trig<cycle_times[:,0])[0][0]
        #cycle_before.extend([c0-2, c0-1])
        #cycle_after.extend([c0, c0+1])
        cycle_before.extend([c0-1,])
        cycle_after.extend([c0, ])
    cycle_before, cycle_after = np.array(cycle_before), np.array(cycle_after)
    cycle_before =cycle_before[(cycle_before>=first_cycle) & (cycle_before<nb_cycles)]
    cycle_after =cycle_after[(cycle_after>=first_cycle) & (cycle_after<nb_cycles)]
    
    cycle_all = np.arange(first_cycle, nb_cycles)
    cycle_other = np.setdiff1d(cycle_all, cycle_before)
    cycle_other = np.setdiff1d(cycle_other, cycle_after)
    
    
    # clean artefact
    cycle_artefact =[]
    for c in cycle_all:
        out = is_artifact_in_window(anasig, cycle_times[c, 0], cycle_times[c+1, 0], dbinfo = dbinfo)
        if out:
            cycle_artefact.append(c)
    cycle_artefact = np.array(cycle_artefact)
    
    cycle_all = np.setdiff1d(cycle_all, cycle_artefact)
    cycle_before = np.setdiff1d(cycle_before, cycle_artefact)
    cycle_after = np.setdiff1d(cycle_after, cycle_artefact)
    cycle_other = np.setdiff1d(cycle_other, cycle_artefact)
    
    return cycle_all, cycle_before, cycle_after, cycle_other
        

def compute_mean_cycle_maps(runs, rc):
    """
    This compute sum of cycles for one electrodes:
       * all, before, after and other
    """
    n1, n2, n3, n4 = 0, 0, 0, 0
    sum_cycle_all, sum_cycle_before, sum_cycle_after, sum_cycle_other = None, None, None, None
    for run in runs:
        seg = run.segments[0]
        respsig = seg.respirationsignals[0]
        anasig = seg.analogsignals.filter_by(channel_index = rc.index).first()
        assert anasig is not None, 'No anasig'
        
        tfr, first_cycle,  nb_cycles, cycle_map, time_to_cycles, times, cycle_times, all_cycles = compute_full_cycle_map(anasig.id, respsig.id, f_start, f_stop, deltafreq, f0, nb_point_by_cycle)
        
        cycle_all, cycle_before, cycle_after, cycle_other = get_cycle_nums(run, anasig, first_cycle, nb_cycles, cycle_times)
        
        m1 = sum_cycles(cycle_map, cycle_all, first_cycle)
        n1 += len(cycle_all)
        if sum_cycle_all is None:
            sum_cycle_all = m1
        else:
            sum_cycle_all += m1
        
        m2 = sum_cycles(cycle_map, cycle_before, first_cycle)
        n2 += len(cycle_before)
        if sum_cycle_before is None:
            sum_cycle_before = m2
        else:
            sum_cycle_before += m2
        
        m3 = sum_cycles(cycle_map, cycle_after, first_cycle)
        n3 += len(cycle_after)
        if sum_cycle_after is None:
            sum_cycle_after = m3
        else:
            sum_cycle_after += m3

        m4 = sum_cycles(cycle_map, cycle_other, first_cycle)
        n4 += len(cycle_other)
        if sum_cycle_other is None:
            sum_cycle_other = m4
        else:
            sum_cycle_other += m4
    
    return { 'cycle_all' : {'sum' : sum_cycle_all, 'count' : n1},
                      'cycle_before': {'sum' : sum_cycle_before, 'count' : n2},
                      'cycle_after': {'sum' : sum_cycle_after, 'count' : n3},
                      'cycle_other': {'sum' : sum_cycle_other, 'count' : n4},
                      'full_map' : (tfr, first_cycle,  nb_cycles, cycle_map, time_to_cycles, times, cycle_times, all_cycles),
                      
                    }
    

def plot_summary(subject_name, run_exp, run_index, rc_name):
    subject = session.query(Subject).filter_by(name = subject_name)[0]
    run = subject.runs.filter_by(exp = run_exp).filter_by(index = run_index)[0]
    runs = [run]
    rc = subject.blocks[0].recordingchannelgroups[0].recordingchannels.filter_by(name = rc_name).first()
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    
    results = compute_mean_cycle_maps(runs, rc)
    tfr, first_cycle,  nb_cycles,  cycle_map, time_to_cycles, times, cycle_times, all_cycles = results['full_map']

    trig_times = [trial.triggered_odor_time for trial in run.trials.order_by(text('`index`'))]
    trig_cycles = []
    for t in trig_times:
        ind  = np.argmin(np.abs(t-times))
        trig_cycles.append(time_to_cycles[ind])

    
    # Plot TFR
    map = np.abs(tfr.map).astype('float32')
    tfr_times =  tfr.times.rescale('s').magnitude
    fig, ax = pyplot.subplots()
    ax.set_title('tfr time')
    plot_map(map, ax, tfr_times, clim = clim, colorbar = True)
    for insp, expi in cycle_times:
        ax.axvline(insp, ls = '-', color = 'w')
        ax.axvline(expi, ls = '--', color = 'w')
    for t in trig_times:
        ax.axvline(t, ls = '-', color = 'r')

    ax.set_title('tfr time')

    # Plot time to cycle
    fig, ax = pyplot.subplots()
    for c in range(nb_cycles):
        ax.axvline(cycle_times[c,0], ls = '-', color = 'c')
        ax.axhline(c, ls = '-', color = 'c')
        ax.axvline(cycle_times[c,1], ls = '--', color = 'c')
        ax.axhline(c+inspi_ratio, ls = '--', color = 'c')
    ax.plot(times, time_to_cycles, lw = 2)


    # Plot cycle map
    fig, ax = pyplot.subplots()
    plot_map(cycle_map, ax, all_cycles, clim = clim, colorbar = True)
    for c in range(nb_cycles):
        ax.axvline(c, ls = '-', color = 'w')
        ax.axvline(c+inspi_ratio, ls = '--', color = 'w')
    for tc in trig_cycles:
        ax.axvline(tc, ls = '-', color = 'r')
        
    
    # Plot  mean cycle map
    fig, axs = pyplot.subplots(ncols = 4)
    axsD = {'cycle_all' : axs[0], 'cycle_before' : axs[1], 'cycle_after' : axs[2], 'cycle_other' : axs[3],}
    for k, ax in axsD.items():
        result = results[k]
        if result['count']==0: continue
        
        mean_cycle_map = result['sum']/result['count']
        plot_map(mean_cycle_map, ax, one_cycle, clim = clim, colorbar = (k=='cycle_all'))
        ax.set_title('{} N={}'.format(k, result['count']))
        ax.axvline(inspi_ratio, ls = '--', color = 'w')
    fig.text(.5, 0.02, '{} {} {}'.format(subject.name, run_names, rc.name))
    
    # cycle envelop
    fig, ax = pyplot.subplots()
    mask = (tfr.freqs.magnitude>f_start_env) & (tfr.freqs.magnitude<=f_stop_env)
    band_freqs = tfr.freqs.magnitude[mask]
    envelop = np.max(cycle_map[:, mask], axis=1)
    #~ envelop_freqs = band_freqs[np.argmax(cycle_map[:, mask], axis=1)]
    ax.plot(all_cycles, envelop, color = 'b')
    for c in range(nb_cycles):
        ax.axvline(c, ls = '-', color = 'c')
        ax.axvline(c+inspi_ratio, ls = '--', color = 'c')
    for tc in trig_cycles:
        ax.axvline(tc, ls = '-', color = 'r')
    ax.set_title('envelop {}-{}Hz'.format(f_start_env, f_stop_env))

    
    #plot by trial envelop
    plot_by_trial_envelop_cycle([run], rc_name)



def plot_and_savefig_mean_cyle_maps(run_ids):
    """
    This compute cycle maps for all electrodes given a run_ids list.
    All figure are save in png.
    This can be lanuch with multiprocessing. (because new open_db at each call)
    
    """
    dbinfo = open_db(url, myglobals = locals(), use_global_session = False, 
                            object_number_in_cache = 0,  numpy_storage_engine = 'sqltable', compress = compress,
                            relationship_lazy = 'dynamic', predefined_classes = classes)
    session = dbinfo.Session()
    #~ subject = session.query(Subject).get(subject_id)
    runs = [session.query(Run).get(run_id) for run_id in run_ids]
    subject = runs[0].subject
    
    base_dirname = os.path.join(mean_cycles_maps_save_dirname, runs[0].subject.name)
    if not os.path.exists(base_dirname):
        os.mkdir(base_dirname)
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    dirname = os.path.join(base_dirname, run_names)
    #~ print dirname
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    
    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels:
        print subject.name, run_names, rc.name
        try:
            fig, axs = pyplot.subplots(ncols = 4, figsize = (16, 12))
            axsD = {'cycle_all' : axs[0], 'cycle_before' : axs[1], 'cycle_after' : axs[2], 'cycle_other' : axs[3],}
            
            results = compute_mean_cycle_maps(runs, rc)
            
            for k, ax in axsD.items():
                result = results[k]
                if result['count']==0: continue
                
                mean_cycle_map = result['sum']/result['count']
                plot_map(mean_cycle_map, ax, one_cycle, clim = clim, colorbar = (k=='cycle_all'))
                ax.set_title('{} N={}'.format(k, result['count']))
                ax.axvline(inspi_ratio, ls = '--', color = 'w')
            
            fig.text(.5, 0.02, '{} {} {}'.format(subject.name, run_names, rc.name))
            fig.savefig(os.path.join(dirname, '{}.png'.format(rc.name)))
            pyplot.close(fig)
        except:
            print '!!!WTF'
    session.close()
    dbinfo.Session.get_bind().dispose()


def plot_and_savefig_mean_cyle_maps_all(with_multiprocess = False, pool_size = 3):
    args_list = []
    for subject in session.query(Subject):
        runs = subject.runs.filter_by(exp = 'E').filter_by(index = 1)[:1]
        args_list.append([run.id for run in runs])
        runs = subject.runs.filter_by(exp = 'E').filter_by(index = 2)[:1]
        args_list.append([run.id for run in runs])
        runs = subject.runs.filter_by(exp = 'R').order_by('`index`')[:3]
        args_list.append([run.id for run in runs])
    if with_multiprocess:
        pool = mp.Pool(pool_size)
        pool.map(plot_and_savefig_mean_cyle_maps, args_list)
        pool.join()
    else:
        for args in args_list:
            plot_and_savefig_mean_cyle_maps(args)

def plot_and_savefig_diff_cycle_maps(run_ids):
    """
    This compute cycle maps for all electrodes given a run_ids list.
    All figure are save in png.
    This can be lanuch with multiprocessing. (because new open_db at each call)
    
    """
    dbinfo = open_db(url, myglobals = locals(), use_global_session = False, 
                            object_number_in_cache = 0,  numpy_storage_engine = 'sqltable', compress = compress,
                            relationship_lazy = 'dynamic', predefined_classes = classes)
    session = dbinfo.Session()
    #~ subject = session.query(Subject).get(subject_id)
    runs = [session.query(Run).get(run_id) for run_id in run_ids]
    subject = runs[0].subject
    
    base_dirname = os.path.join(diff_cycle_maps_save_dirname, runs[0].subject.name)
    if not os.path.exists(base_dirname):
        os.mkdir(base_dirname)
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    dirname = os.path.join(base_dirname, run_names)
    #~ print dirname
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels.filter_by(name = 'j2'):
        print subject.name, run_names, rc.name
        if 1:
        #~ try:
            results = compute_mean_cycle_maps(runs, rc)
            if results['cycle_before']['count']==0: continue
            elif results['cycle_other']['count']==0: continue
            elif results['cycle_after']['count']==0: continue
            fig, axs = pyplot.subplots(ncols = 2, figsize = (16, 12))
            axsD = {'cycle_before_other' : axs[0], 'cycle_after_other' : axs[1],}
            Diff_TF_all_cycles_before = results['cycle_before']['sum']/results['cycle_before']['count']\
                    - results['cycle_other']['sum']/results['cycle_other']['count']
            plot_map(Diff_TF_all_cycles_before, axs[0], one_cycle, clim = [0, 5], colorbar = True)
            axs[0].set_title('cycle_before_other')
            axs[0].axvline(inspi_ratio, ls = '--', color = 'w')
            
            Diff_TF_all_cycles_after = results['cycle_after']['sum']/results['cycle_after']['count']\
                                    - results['cycle_other']['sum']/results['cycle_other']['count']
            
            plot_map(Diff_TF_all_cycles_after, axs[1], one_cycle, clim = [0, 5], colorbar = False)
            axs[1].set_title('cycle_after_other')
            axs[1].axvline(inspi_ratio, ls = '--', color = 'w')

            fig.text(.5, 0.02, '{} {} {}'.format(subject.name, run_names, rc.name))
            
            fig.savefig(os.path.join(dirname, '{}.png'.format(rc.name)))
            pyplot.close(fig)
        #~ except:
            #~ print '!!!WTF'
    session.close()
    dbinfo.Session.get_bind().dispose()


def plot_and_savefig_diff_cycle_maps_all(with_multiprocess = False, pool_size = 2):
    args_list = []
    for subject in session.query(Subject).filter_by(name = 'LEFC'):
        #runs = subject.runs.filter_by(exp = 'E').filter_by(index = 1)[:1]
        runs = subject.runs.filter_by(exp = 'E').order_by('index')[:2]
        #~ runs = subject.runs.filter_by(exp = 'E').filter_by(index = 2)[:1]
        #~ args_list.append([run.id for run in runs])
        #runs = subject.runs.filter_by(exp = 'R').order_by('`index`')[:3]
        args_list.append([run.id for run in runs])
    if with_multiprocess:
        pool = mp.Pool(pool_size)
        pool.map(plot_and_savefig_diff_cycle_maps, args_list)
        pool.join()
    else:
        for args in args_list:
            plot_and_savefig_diff_cycle_maps(args)


def  plot_by_trial_envelop_cycle(runs, rc_name):
    
    #runs = [session.query(Run).get(run_id) for run_id in run_ids]
    
    all_env = []
    all_freqs = []
    for run in runs:
        rc = run.subject.blocks[0].recordingchannelgroups[0].recordingchannels.filter_by(name = rc_name).first()
        seg = run.segments[0]
        respsig = seg.respirationsignals[0]
        anasig = seg.analogsignals.filter_by(channel_index = rc.index).first()
        tfr, first_cycle,  nb_cycles, cycle_map, time_to_cycles, times, cycle_times, all_cycles = compute_full_cycle_map(anasig.id, respsig.id, f_start, f_stop, deltafreq, f0, nb_point_by_cycle)
        
        cycle_all, cycle_before, cycle_after, cycle_other = get_cycle_nums(run, anasig, first_cycle, nb_cycles, cycle_times)
        cycle_before = cycle_before[np.in1d(cycle_before+1, cycle_after)]
        
        mask = (tfr.freqs.magnitude>f_start_env) & (tfr.freqs.magnitude<=f_stop_env)
        band_freqs = tfr.freqs.magnitude[mask]
        envelop = np.max(cycle_map[:, mask], axis=1)
        envelop_freqs = band_freqs[np.argmax(cycle_map[:, mask], axis=1)]
        
        for c in cycle_before:
            c2 = c - first_cycle
            env = envelop[c2*nb_point_by_cycle:(c2+2)*nb_point_by_cycle]
            all_env.append(env)
            
            all_freqs.append(envelop_freqs[c2*nb_point_by_cycle:(c2+2)*nb_point_by_cycle])
            
            
    all_env = np.array(all_env)
    all_freqs = np.array(all_freqs)
    
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    
    fig, axs = pyplot.subplots(ncols = 3)
    
    axs[0].plot(two_cycle, all_env.transpose(), color = 'k')
    axs[0].set_ylim(clim[0], clim[1]*2)
    axs[0].set_ylabel('envelop amplitude')
    axs[0].set_title('{} envelop {}-{}Hz'.format(rc_name, f_start_env, f_stop_env))

    axs[1].plot(two_cycle, all_freqs.transpose(), color = 'k')
    axs[1].set_ylim(clim)
    axs[1].set_ylabel('envelop frequency')
    
    
    im = axs[2].imshow(all_env,
                            interpolation='nearest', 
                            extent=(0, 2, 0, cycle_before.size),
                            origin ='lower' ,
                            aspect = 'auto')
    im.set_clim(clim)
    fig.colorbar(im)
    
    
    
    for ax in axs:
        ax.axvline(1., ls = '-', color = 'c')
        ax.axvline(inspi_ratio, ls = '--', color = 'c')
        ax.axvline(inspi_ratio+1., ls = '--', color = 'c')
        ax.set_xlabel('cycle before and after')
    
    subject = runs[0].subject
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    fig.text(.5, 0.02, '{} {} {}'.format(subject.name, run_names, rc.name))        
    
    return fig


def plot_and_savefig_by_trial_envelop_cycle(run_ids):
    dbinfo = open_db(url, myglobals = locals(), use_global_session = False, 
                            object_number_in_cache = 0,  numpy_storage_engine = 'sqltable', compress = compress,
                            relationship_lazy = 'dynamic', predefined_classes = classes)
    session = dbinfo.Session()
    runs = [session.query(Run).get(run_id) for run_id in run_ids]
    subject = runs[0].subject


    base_dirname = os.path.join(bytrial_envelop_cycles_save_dirname, runs[0].subject.name)
    if not os.path.exists(base_dirname):
        os.mkdir(base_dirname)
    run_names =  ' '.join(['{}{}'.format(run.exp, run.index) for run in runs])
    dirname = os.path.join(base_dirname, run_names)
    #~ print dirname
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels.filter_by():
        print subject.name, run_names, rc.name
        try:
        #~ if 1:
            fig = plot_by_trial_envelop_cycle(runs, rc.name)
            fig.set_figsize = (16, 12)
            fig.savefig(os.path.join(dirname, '{}.png'.format(rc.name)))
            pyplot.close(fig)
        except:
            print '!!!WTF'
    session.close()    
    dbinfo.Session.get_bind().dispose()


def plot_and_savefig_by_trial_envelop_cycle_all(with_multiprocess = False, pool_size = 3):
    args_list = []
    for subject in session.query(Subject):
        runs = subject.runs.filter_by(exp = 'E')[:2]
        args_list.append([run.id for run in runs])
        #~ runs = subject.runs.filter_by(exp = 'E').filter_by(index = 2)[:1]
        #~ args_list.append([run.id for run in runs])
        #~ runs = subject.runs.filter_by(exp = 'R').order_by('`index`')[:3]
        #~ args_list.append([run.id for run in runs])
    if with_multiprocess:
        pool = mp.Pool(pool_size)
        pool.map(plot_and_savefig_by_trial_envelop_cycle, args_list)
        pool.join()
    else:
        for args in args_list:
            plot_and_savefig_by_trial_envelop_cycle(args)

def grand_mean_cycle_maps():
    sheetname = 'mHipp'
    
    conditions = ['cycle_all', 'cycle_before', 'cycle_after', 'cycle_other']
    
    
    mean_list = pd.read_excel('mean_elec_cycle_maps.xlsx', sheetname=sheetname, index_col = None)
    all_mean = { k : [] for k in  conditions }
    for i, row in mean_list.iterrows():
        #~ print row
        subject = session.query(Subject).filter_by(name =  row['subject_name']).first()
        runs = subject.runs.filter_by(exp = row['run_exp']).filter_by(index = row['run_index']).all()
        rc = subject.blocks[0].recordingchannelgroups[0].recordingchannels.filter_by(name = row['rc_name']).first()
        if rc is None:
            print row['subject_name'], row['rc_name']
            continue
        results = compute_mean_cycle_maps(runs, rc)
        
        for k in conditions:
            result = results[k]
            mean_cycle_map = result['sum']/result['count']
            all_mean[k].append(mean_cycle_map)
    
    fig, axs = pyplot.subplots(ncols = len(conditions), figsize = (16, 12))    
    axsD = { cond : axs[i] for i,cond in enumerate(conditions) }    
    for k in conditions:
        all = np.array(all_mean[k])
        m = np.mean(all, axis=0)
        ax = axsD[k]
        plot_map(m, ax, one_cycle, clim = clim, colorbar = (k=='cycle_all'))
        
        ax.set_title(u'{} {} N={}'.format(sheetname, k, all.shape[0]))
        ax.axvline(inspi_ratio, ls = '--', color = 'w')
    
    fig.savefig(sheetname+'.png')
    
    
    
    
if __name__ == '__main__':
    #estimate_inspi_ratio()
    #plot_summary('PIRJ', 'E', '2', "t3")
    #~ pyplot.show()
    
    #~ plot_and_savefig_mean_cyle_maps_all(with_multiprocess = False)
    #~ plot_and_savefig_mean_cyle_maps_all(with_multiprocess = True, pool_size = 6)
    #~ mean_cyle_maps_cule_maps()
    
    #plot_and_savefig_diff_cycle_maps_all(with_multiprocess = False)
    #plot_and_savefig_diff_cycle_maps_all(with_multiprocess = True, pool_size = 4)
    
    plot_by_trial_envelop_cycle([session.query(Run).get(354)], "d'11")
    pyplot.show()
    
    #plot_and_savefig_by_trial_envelop_cycle_all(with_multiprocess = False)
    #~ plot_and_savefig_by_trial_envelop_cycle_all(with_multiprocess = True, pool_size = 3)
    
    #~ grand_mean_cycle_maps()
    #~ pyplot.show()
    
    
