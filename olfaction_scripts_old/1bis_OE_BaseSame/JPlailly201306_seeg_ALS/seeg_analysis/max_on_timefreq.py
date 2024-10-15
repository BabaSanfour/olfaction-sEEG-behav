# -*- coding: utf-8 -*-


import sys, os
sys.path.append('../behavior') 

from connection import *

import params

from OpenElectrophy.timefrequency import TimeFreq


from artifact_detection import is_artifact_in_window,clean_trig

import scipy.signal

import quantities as pq
import numpy as np 

from matplotlib import pyplot
from matplotlib.cm import get_cmap
from joblib import Memory

memory  = Memory(cachedir=params.joblib_cachedir)


#~ from tools import get_envelope_in_band, get_timefreq
from timefreqtools import*



def max_on_timefreq():
    subject = session.query(Subject).filter_by(name = 'LEFC').first()
    #~ run = subject.runs.filter_by(exp= 'E').filter_by(index=2).first()
    run = subject.runs.filter_by(exp= 'E').first()
    print run.subject
    print run
    rcg = subject.blocks[0].recordingchannelgroups.filter_by(name = 'group a').first()    
    assert rcg is not None, 'pas de group'
    
    rcs = rcg.recordingchannels.all()
    rcs = rcs[:3]
    print rcs
    #~ rcs = rcs[:2]
    
    n = len(rcs)
    
    fig, ax1s = pyplot.subplots(nrows=n, ncols=3, figsize=[16,8])
    fig.text(.5, .05, run.filename)
    
    trig_times = [ trial.first_inspiration_time for trial in run.trials]
    #~ trig_times = [ trial.triggered_odor_time for trial in run.trials]
    

    for i, rc in enumerate(rcs):
        print i, 'Channel', rc.index, rc.name
        print rc
        
        
        #find the signal
        anasig = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc.id)[0]
        print anasig
        neosig = anasig.to_neo()
        
        t1 = -1.
        t2 = 3.

        f_start = 1.
        f_stop = 35.

        #~ f_start = params.f_start
        #~ f_stop = params.f_stop
        
        #~ tfr = get_timefreq(anasig.id, f_start, f_stop, (f_stop-f_start)/40)
        tfr = computations.get_timefreq(anasig.id, f_start, f_stop, f0 = 2.5, deltafreq = (f_stop-f_start)/40)
        map = np.abs(tfr.map)
        times = tfr.times.rescale('s').magnitude
        freqs = tfr.freqs.rescale('Hz').magnitude
        
        all_max = [ ]
        for j, trig_time in enumerate(trig_times):
            if is_artifact_in_window(anasig, trig_time+t1, trig_time+t2):
                print 'SKIP'
                continue
            
            sel = (times>=trig_time+t1) & (times<trig_time+t2)
            small = map[sel, :]
            ind_t, ind_f =  np.unravel_index(small.argmax(), map.shape)
            
            ampl_max = small[ind_t, ind_f]
            freq_max = freqs[ind_f]
            offset_max = times[sel][ind_t] - trig_time
            #~ print ampl_max, freq_max, offset_max
            all_max.append([ampl_max, freq_max, offset_max])
        all_max = np.array(all_max)
        
        
        labels = ['amplitude_max', 'freq_max', 'offset_max']
        
        for j in range(3):
            if j==0:
                bins = np.arange(0,50,2)
            elif j==1:
                bins = np.arange(f_start,f_stop,(f_stop-f_start)/40)
            elif j==2:
                bins = np.arange(t1,t2, (t2-t1)/40)
                
            ax = ax1s[i,j]
            y, x = np.histogram(all_max[:,j], bins= bins)
            ax.plot(x[:-1], y)
            
            if j!=0: ax.set_yticks([])
            if i!=len(rcs)-1: ax.set_xticks([])
            if i==0 :
                ax.set_title(labels[j])
            if j==0 :
                ax.set_ylabel(rc.name)
        
    pyplot.show()
            
    
    
    
    
    
    
    
if __name__ == '__main__':
    max_on_timefreq()
    