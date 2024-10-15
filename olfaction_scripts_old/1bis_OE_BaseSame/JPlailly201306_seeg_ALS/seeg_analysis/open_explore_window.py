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
memory  = Memory(cachedir=params.joblib_cachedir)
memory  = Memory(params.joblib_cachedir, mmap_mode = 'r')
timefreqtools.cache_all_computation(memory)




def get_data():
    subject = session.query(Subject).filter_by(name = 'LEFC').first()

    trials = [ ]
    for run in  subject.runs.order_by('exp').order_by('`index`'):
        for trial in run.trials:
            if trial.run.segments.first() is None: continue
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
    return subject,  trials, rcgs
    


def open_explore_window():
    subject, trials, rcgs = get_data()
    trial_sub_list = [{ 'label' : 'All trials', 'trials' : trials}]
    app = pg.mkQApp()
    w = ExploreTimeFreqWindow(title = u'Explore freqs for {}'.format(subject.name),
                    trial_sub_list = trial_sub_list,
                    rcgs = rcgs,
                    ylim = 50,
                    )
    w.show()
    app.exec_()



def test1():
    subject, trials, rcgs = get_data()
    trials = trials[:3]
    rcs = subject.blocks[0].recordingchannelgroups.filter_by(name = 'Group b').first().recordingchannels
    
    rcs = rcs[3:4]
    
    print rcs[0].analogsignals[0].sampling_rate
    
    print dir(timefreqtools)
    #~ timefreqtools.PlotSignalOnTrig(rcs = rcs, trials = trials, t_start = -1., t_stop = 2., trig_alignement = 'triggered_odor_time')
    #~ timefreqtools.PlotFilteredSignalOnTrig(rcs = rcs, trials = trials, t_start = -1., t_stop = 2., f_start = 15., f_stop=35., envelop_method = 'timefreq',)
    #~ timefreqtools.PlotTimeFreqOnTrig(rcs = rcs, trials = trials, t_start = -1., t_stop = 2., f_start = 15., f_stop=35., nb_line = 40., ylim =50)
    timefreqtools.PlotMeanTimeFreqOnTrig(rcs = rcs, trials = trials, t_start = -1., t_stop = 2., f_start = 15., f_stop=35., nb_line = 40., ylim =50)
    
    
    #~ timefreqtools.PlotEnvelopOnTrig(rcs = rcs, trials = trials, t_start = -1., t_stop = 2., f_start = 15., f_stop=35.,envelop_method = 'timefreq',)
    
    
    #~ timefreqtools.PlotLongTermTimeFreqProbaDensity(rcs = rcs, trials = trials, f_start = 15., f_stop=512./2., ylim = 50, nb_line = 80)
    #~ timefreqtools.PlotWindowedFreqProbaDensity(rcs = rcs, trials = trials,  f_start = 1., f_stop=512./2., ylim = 50, nb_line = 80)
    
    
    timefreqtools.show()

def test2():
    subject, trials, rcgs = get_data()
    trials = trials[:3]
    rcs = subject.blocks[0].recordingchannelgroups.filter_by(name = 'Group b').first().recordingchannels[:2]

    app = pg.mkQApp()
    w = MultiFigure(rcs = rcs, trials = trials,
                    t_start = -1., t_stop = 2.,
                    f_start = 15., f_stop=35.,
                    nb_line = 40.,
                    envelop_method = 'timefreq',
                    
                    
                    #~ plot_selection = ['PlotSignalOnTrig', 'PlotTimeFreqOnTrig', ],
                    )
    w.show()
    
    app.exec_()


if __name__ == '__main__':
    #~ test1()
    #~ test2()
    
    open_explore_window()



