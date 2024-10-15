# -*- coding: utf-8 -*-


import sys, os
sys.path.append('../behavior') 
from connection import *
import params


import pyqtgraph as pg
from timefreqtools import ExploreTimeFreqWindow, MultiFigure
import timefreqtools


from artifact_detection import is_artifact_in_window

timefreqtools.set_session(session, dbinfo)

from joblib import Memory
memory  = Memory(cachedir=params.joblib_cachedir)
memory  = Memory(params.joblib_cachedir, mmap_mode = 'r')
timefreqtools.cache_all_computation(memory)



t_start = -1.
t_stop = 2.
f_start = 15.
f_stop=35.
nb_line = 40.
ylim =50


def one_tfr_mean(subject_name = 'LEFC', rc_index = 9):
    subject = session.query(Subject).filter_by(name = subject_name).first()

    rcg_all = subject.blocks[0].recordingchannelgroups[0]
    rc = rcg_all.recordingchannels.filter_by(index = rc_index).first()
    
    assert rc is not None, 'pas de rc'
    rcs = [rc]


    trials = [ ]
    for run in  subject.runs.order_by('exp').order_by('`index`'):
        for trial in run.trials:
            seg = trial.run.segments.first()
            if seg is None: continue
            
            mon_trig = trial.triggered_odor_time
            #~ mon_trig = trial.time
            #~ mon_trig = trial.recognition_time
            
            
            #remove artefact
            anasig = seg.analogsignals.filter_by(recordingchannel_id = rc.id).one()
            if is_artifact_in_window(anasig, mon_trig+t_start, mon_trig+t_stop):
                print 'artefact on ', subject_name, rc_index, trial.index
                continue
            
            trial._info  = { 'label':u'{}{} trial{}  {} {}'.format(trial.run.exp, trial.run.index, 
                                                trial.index, trial.score_recognition, trial.score_episodic_strict),
                                        'events' : [('mon_trig', mon_trig, 'r'), ],
                                        'epochs' : [ ],
                                        'segment' : trial.run.segments.first(),
                                    }
            trials.append(trial)
    
    analyse_plot = timefreqtools.PlotMeanTimeFreqOnTrig(rcs = rcs, trials = trials, 
                        t_start = t_start, t_stop = t_stop, 
                        f_start = f_start, f_stop=f_stop,
                        nb_line = nb_line, ylim =ylim,
                        trig_alignement = 'mon_trig')
    
    
    analyse_plot.fig.savefig('mean_time_maps/{} {}.png'.format(subject.name, rc.name))
    

def all_tfr_mean():
    for subject in session.query(Subject):
        
        rcg_all = subject.blocks[0].recordingchannelgroups[0]
        
        for rc in rcg_all.recordingchannels:
            print 'Compute', subject.name, rc.index
            
            try:
                one_tfr_mean(subject_name = subject.name, rc_index = rc.index)
            except:
                print '    ##### ERROR', subject.name, rc.index





if __name__ =='__main__':
    #one_tfr_mean()
    
    all_tfr_mean()
    
    