# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 
from connection import *

import quantities as pq
import numpy as np
import scipy.signal
from OpenElectrophy.timefrequency import TimeFreq


from matplotlib import pyplot
import time
import neo

import scipy.interpolate


f_start = 1.
f_stop = 100.
deltafreq = .5
f0 = 2.5
sampling_rate_tfr = f_stop*3



nb_point_by_cycle = 30
inspi_ratio = .4
nb_point_inspi = int(nb_point_by_cycle*inspi_ratio)
nb_point_expi = nb_point_by_cycle - nb_point_inspi



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


def resp_cycle_frequency_map():
    subject = session.query(Subject).filter_by(name = 'LEFC')[0]
    run = subject.runs.order_by('date')[0]
    print run
    seg = run.segments[0]
    anasig = seg.analogsignals.filter_by(name = 'M11')[0]
    respsig = seg.respirationsignals[0]
    cycle_times = respsig.cycle_times.magnitude
    print cycle_times.shape
    nb_cycles = cycle_times.shape[0]-1


    
    
    neosig = anasig.to_neo()
    
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
    print map.shape
    times =  tfr.times.rescale('s').magnitude
    
    
    fig, ax = pyplot.subplots()
    ax.set_title('tfr time')
    #~ tfr.plot(ax = ax)
    plot_map(map, ax, times, clim = None, colorbar = True)
    
    for insp, expi in cycle_times:
        ax.axvline(insp, ls = '-', color = 'w')
        ax.axvline(expi, ls = '--', color = 'w')
    
    keep = (times>=cycle_times[0,0]) & (times<cycle_times[-1,0])
    map2 = map[keep, :]
    times2 = times[keep]

    fig, ax = pyplot.subplots()
    ax.set_title('tfr time')
    #~ tfr.plot(ax = ax)
    plot_map(map2, ax, times2, clim = None, colorbar = True)
    
    for insp, expi in cycle_times:
        ax.axvline(insp, ls = '-', color = 'w')
        ax.axvline(expi, ls = '--', color = 'w')
    
    # construct cycle_step
    cycle_steps = np.zeros(times2.shape)*np.nan
    for c in range(nb_cycles):
        #~ #One segment
        #~ mask = (times2>=cycle_times[c, 0]) & (times2<cycle_times[c+1, 0])
        #~ cycle_steps[mask] = np.linspace(c, c+1, num = np.sum(mask))
        
        #2 segments : inspi + expi
        mask = (times2>=cycle_times[c, 0]) & (times2<cycle_times[c, 1])
        cycle_steps[mask] = np.linspace(c, c+inspi_ratio, num = np.sum(mask))
        mask = (times2>=cycle_times[c, 1]) & (times2<cycle_times[c+1, 0])
        cycle_steps[mask] = np.linspace(c+inspi_ratio, c+1, num = np.sum(mask))
    
    fig, ax = pyplot.subplots()
    for c in range(nb_cycles):
        ax.axvline(cycle_times[c,0], ls = '-', color = 'c')
        ax.axhline(c, ls = '-', color = 'c')
        ax.axvline(cycle_times[c,1], ls = '--', color = 'c')
        ax.axhline(c+inspi_ratio, ls = '--', color = 'c')
        
    ax.plot(times2, cycle_steps, lw = 2)
    
    
    interp = scipy.interpolate.interp1d(cycle_steps, map2, kind='linear', axis=0)
    cycles = np.arange(0, nb_cycles, 1./nb_point_by_cycle)
    print cycle_steps[0], cycle_steps[-1]
    print cycles[0], cycles[-1]

    cycle_map = interp(cycles)
    print cycle_map.shape
    
    
    fig, ax = pyplot.subplots()
    plot_map(cycle_map, ax, cycles, clim = None, colorbar = True)
    for c in range(nb_cycles):
        ax.axvline(c, ls = '-', color = 'w')
        ax.axvline(c+inspi_ratio, ls = '--', color = 'w')
    
    all_cycle_map = np.zeros((nb_point_by_cycle, map2.shape[1]))
    for c in range(nb_cycles):
        all_cycle_map += cycle_map[c*nb_point_by_cycle:(c+1)*nb_point_by_cycle, :]
    all_cycle_map /= nb_cycles
    fig, ax = pyplot.subplots()
    plot_map(all_cycle_map, ax, np.linspace(0,1,nb_point_by_cycle), clim = None, colorbar = True)
    
    
    
        
    
    
    

    
if __name__ == '__main__':
    resp_cycle_frequency_map()
    pyplot.show()
