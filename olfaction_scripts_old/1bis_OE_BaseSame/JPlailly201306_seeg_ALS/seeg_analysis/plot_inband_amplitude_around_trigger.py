# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 
from connection import *
import params


import pyqtgraph as pg
from timefreqtools import ExploreTimeFreqWindow, MultiFigure
import timefreqtools


timefreqtools.set_session(session, dbinfo)

from joblib import Memory
#~ memory  = Memory(cachedir=params.joblib_cachedir)
memory  = Memory(params.joblib_cachedir, mmap_mode = 'r')
timefreqtools.cache_all_computation(memory)



params = dict(
            t_start = -1.5,
            t_stop = 3.,
            f_start = 15.,
            f_stop = 40.,
            baseline_t1 = -1.5,
            baseline_t2 = -0.5,

            #~ baseline_mode = 'no-normalisation'
            #~ baseline_mode = 'ratio'
            #~ baseline_mode = 'difference'
            #~ baseline_mode = 'z-score'

            #~ baseline_strategy = 'by-trial'
            #~ baseline_strategy = 'by-run'
            
            ylim = 100,

            target_t1 = 2.,
            target_t2 = 3.,
            nb_line = 40.,
            f0 = 2.5,
            n_max = 3,
            envelop_method = 'timefreq',
)



def plot_amplitude_around_trigger():
    subject = session.query(Subject).filter_by(name = 'SEMC').first()

    trials = [ ]
    for run in  subject.runs.order_by('exp').order_by('`index`'):
        for trial in run.trials:
            trial._info  = { 'label':u'{}{} trial{}  {} {}'.format(trial.run.exp, trial.run.index, 
                                                trial.index, trial.score_recognition, trial.score_episodic_strict),
                                        #~ 'time' : pipet.triggered_time,
                                        'events' : [('triggered_odor_time', trial.triggered_odor_time, 'r'), 
                                                                ('on_first_inspi', trial.first_inspiration_time, '#00BFFF'),
                                                                #~ ('wanted_odor_time', trial.wanted_odor_time, 'r'),
                                                                ('on_trial_time', trial.time, 'w'),
                                                                ('recognition_time', trial.recognition_time, 'm'),
                                                                ('context_time', trial.context_time, 'g'),
                                                                
                                                            ],
                                        'epochs' : [ ],
                                        'segment' : trial.run.segments.first(),
                                    }
            trials.append(trial)
    rcgs = subject.blocks[0].recordingchannelgroups.all()
    rcs = subject.blocks[0].recordingchannelgroups.filter_by(name = 'Group b').first().recordingchannels
    
    trials = trials[15:18]
    rcs = rcs[3:4]



    timefreqtools.PlotSignalOnTrig(rcs = rcs, trials = trials,**params)
    timefreqtools.PlotFilteredSignalOnTrig(rcs = rcs, trials = trials, **params)
    timefreqtools.PlotTimeFreqOnTrig(rcs = rcs, trials = trials, **params)
    timefreqtools.PlotEnvelopOnTrig(rcs = rcs, trials = trials, **params)
    timefreqtools.PlotEnvelopOnTrigSummary(rcs = rcs, trials = trials, **params)
    
    timefreqtools.show()
    
    
    
    

if __name__ == '__main__':
    
    plot_amplitude_around_trigger()
    
    
    
