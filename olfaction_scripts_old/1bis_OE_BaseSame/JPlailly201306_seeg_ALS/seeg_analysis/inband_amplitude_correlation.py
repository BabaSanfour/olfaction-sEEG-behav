# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 

from connection import *

from params import base_dir,joblib_cachedir
from params import subject_ids,trial_R_conds

from params import t_win_start,t_win_stop,baseline_t_win_start,baseline_t_win_stop,f_start,f_stop,envelop_method,baseline_mode


from OpenElectrophy.timefrequency import TimeFreq

from inband_amplitude_around_trigger import get_envelope_in_band

from artifact_detection import is_artifact_in_window,detect_artifact_on_one_sig,return_artifact_timings,clean_trig

from dmgraphanalysis_nodes.utils_stats import return_signif_code

from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat,plot_int_mat
        
import itertools

#~ from artifact_detection import is_artifact_in_window

import scipy.signal
from scipy.stats import ttest_ind
from matplotlib import pyplot

from joblib import Memory

memory  = Memory(cachedir=joblib_cachedir)

def export_rcs_to_np(run,rcs):
    
    sigs = []
    
    ampl_sigs = []
    
    if len(rcs) != 0:
        
        for i in range(len(rcs)):
            anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
            t, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
            
            neosig_i = anasig_i.to_neo()
            sig_i = neosig_i.magnitude
            
            times = neosig_i.times.rescale('s').magnitude
            
            #print sig_i
            
            sigs.append(sig_i)
            ampl_sigs.append(ampl_i)
            
    else:
        
        times = np.array([])
    #print sigs
    
    return np.array(sigs,dtype = float),np.array(ampl_sigs,dtype = float),times
    
def return_keep_all_non_art(run,rcs, timings):
    
    keep = np.ones(shape = (len(rcs),timings.shape[0]), dtype = 'bool')
    
    #print keep.shape
    
    #### computing artefacts
    for i,rc in enumerate(rcs):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc.id)[0]
        
        art_start,art_stop = return_artifact_timings(anasig_i)
        
        keep[i,:] = compute_keep_from_timings(timings,zip(art_start,art_stop))
        
        #print keep[i,:]
        
    #print keep.shape
    
    return np.array( np.sum(keep,axis = 0) == 0,dtype = 'bool')
    
def compute_keep_from_timings(timings,art_starts_stops):

    keep = np.zeros(shape = timings.shape, dtype = 'bool')
    
    for art_seq in art_starts_stops:
        
        keep[(timings >= art_seq[0]) & (timings < art_seq[1])] = True
        
    return keep
    
def return_correct_non_trig_ampl_sigs(run,rcs):
    
    ##### transfer data to np array
    np_sigs,np_ampl_sigs,timings = export_rcs_to_np(run,rcs)
    
    keep_all_non_art = return_keep_all_non_art(run,rcs,timings)
    
    print np.sum(keep_all_non_art == True),np.sum(keep_all_non_art == False)
    #print np.sum(keep_all == True)/float(keep_all.shape[0])
    
    #### computing from trigs
    trig_times = np.array([ trial.triggered_odor_time for trial in run.trials ],dtype = 'float')
    
    keep_non_trig = np.logical_not(compute_keep_from_timings(timings,zip(trig_times+t_win_start,trig_times+t_win_stop)))
    
    print keep_non_trig
    print np.sum(keep_non_trig == True),np.sum(keep_non_trig == False)
    
    keep_sigs = np_sigs[:,np.logical_and(keep_all_non_art,keep_non_trig)]
    
    print keep_sigs.shape
    
    keep_ampl_sigs= np_ampl_sigs[:,np.logical_and(keep_all_non_art,keep_non_trig)]
    
    return keep_sigs,keep_ampl_sigs
        
def return_correct_ampl_by_trig_duration(write_path,run,rcs,trial_type = 'all'):
    
    ##### transfer data to np array
    np_sigs,np_ampl_sigs,timings = export_rcs_to_np(run,rcs)
    
    keep_all_non_art = return_keep_all_non_art(run,rcs,timings)
    
    print np.sum(keep_all_non_art == True),np.sum(keep_all_non_art == False)
    #print np.sum(keep_all == True)/float(keep_all.shape[0])
    
    if trial_type == 'all':
        
        trig_times = [ trial.triggered_odor_time for trial in run.trials ]
        
        recognition_times = [ trial.time + trial.recognition_time for trial in run.trials]
        
    else:
        odor_times = [ trial.first_inspiration_time for trial in run.trials if trial.score_recognition == trial_type ]
        
        recognition_times = [ trial.time + trial.recognition_time for trial in run.trials if trial.score_recognition == trial_type ]
    
    for i,art_seq in enumerate(zip(odor_times,recognition_times)):
        
        sel = (timings > art_seq[0]) & (timings <= art_seq[1])
        
        if np.sum(keep_all_non_art[sel]) == np.sum(sel == True):
        
            print art_seq[0],art_seq[1]
                
            trial_ampl =  np_sigs[:,sel]
            
            print trial_ampl.shape
            
            trial_ampl_file = os.path.join(write_path,'ts_' + trial_type+'_ampl_by_odor_trial_' + str(i) + '.npy')
            
            np.save(trial_ampl_file,trial_ampl)

def return_correct_ampl_by_trig_duration_and_rest(write_path,run,rcs,trial_type = 'all'):
    
    ##### transfer data to np array
    np_sigs,np_ampl_sigs,timings = export_rcs_to_np(run,rcs)
    
    keep_all_non_art = return_keep_all_non_art(run,rcs,timings)
    
    #print np.sum(keep_all_non_art == True),np.sum(keep_all_non_art == False)
    
    if trial_type == 'all':
        
        trig_times = [ trial.triggered_odor_time for trial in run.trials ]
        
        recognition_times = [ trial.time + trial.recognition_time for trial in run.trials]
        
    else:
        odor_times = [ trial.first_inspiration_time for trial in run.trials if trial.score_recognition == trial_type ]
        
        recognition_times = [ trial.time + trial.recognition_time for trial in run.trials if trial.score_recognition == trial_type ]
    
        start_rest = [ trial.time for trial in run.trials if trial.score_recognition == trial_type ]
        
        stop_rest = [ trial.time+3.0 for trial in run.trials if trial.score_recognition == trial_type ]
        
    count = 0
    
    for i,art_seq in enumerate(zip(odor_times,recognition_times)):
        
        sel_odor = (timings > art_seq[0]) & (timings <= art_seq[1])
        
        sel_rest = (timings > start_rest[i]) & (timings <= stop_rest[i])
        
        if np.sum(keep_all_non_art[sel_odor]) == np.sum(sel_odor == True) and np.sum(keep_all_non_art[sel_rest]) == np.sum(sel_rest == True):
        
            print art_seq[0],art_seq[1]
                
            odor_ampl =  np_sigs[:,sel_odor]
            
            print odor_ampl.shape
            
            odor_ampl_file = os.path.join(write_path,'ts_' + trial_type+'_ampl_by_odor_trial_' + str(count) + '.npy')
            
            np.save(odor_ampl_file,odor_ampl)
            
            rest_ampl = np_sigs[:,sel_rest]
            
            print rest_ampl.shape
            
            rest_ampl_file = os.path.join(write_path,'ts_' + trial_type+'_ampl_by_rest_trial_' + str(count) + '.npy')
            
            np.save(rest_ampl_file,rest_ampl)
            
            count = count + 1
            
            
            
def return_sig_on_start_stop(t, ampl_i, clean_start_times,clean_stop_times, baseline_mode):

    trig_sigs = []
    
    for start_time,stop_time in zip(clean_start_times,clean_stop_times):
        
        print ampl_i.shape
        
        print start_time,stop_time
        
        sel = (t>=start_time) & (t<stop_time)
        
        print np.sum(sel == True)
        
        if baseline_mode == 'no-normalisation':
            local_sig = ampl_i[sel]
            
        else:
            print "&&&&&&&&&&&& Warning, only no-normalisation is defined &&&&&&&&&&&&&&&&&&&&&&&&&&&"
            
        #print local_sig.shape
        trig_sigs.append(local_sig)
        
    return np.array(trig_sigs,dtype = float)
        
        
        
        
        
def return_sig_on_trig(t, sig, trig_times, baseline_mode):
    
    trig_sigs = []
    
    for i, trig_time in enumerate(trig_times):
        
        sel = (t>=trig_time+t_win_start) & (t<trig_time+t_win_stop)
        
        # this is the baseline window
        sel2 = (t>=trig_time+baseline_t_win_start) & (t<trig_time+baseline_t_win_stop)
        baseline_sig = sig[sel2]
        
        if baseline_mode == 'no-normalisation':
            local_sig = sig[sel]
        elif baseline_mode == 'ratio':
            local_sig = sig[sel]/baseline_sig.mean()
        elif baseline_mode == 'difference':
            local_sig = sig[sel]-baseline_sig.mean()
        elif baseline_mode == 'z-score':
            local_sig = (sig[sel]-baseline_sig.mean())/baseline_sig.std()
        
        #print local_sig.shape
        trig_sigs.append(local_sig)
        
    #print trig_sigs
    
    np_trig_sigs = np.array(trig_sigs,dtype = float)
    
    #print np_trig_sigs.shape
    
    return np_trig_sigs
        
        
############################################################################### trig ampl based ##############################################
def return_all_trig_ampl(run,rcs,trial_type = 'all'):
    
    if trial_type == 'all':
        trig_times = [ trial.triggered_odor_time for trial in run.trials ]
        
    else:
        trig_times = [ trial.triggered_odor_time for trial in run.trials if trial.score_recognition == trial_type ]
    
    
    #print trig_times
    print len(trig_times)
    
    c_trig_times = clean_trig_all_channels(trig_times,run,rcs)
    
    #print c_trig_times
    print len(c_trig_times)
    
    all_trig_ampl = []
    
    
    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        t, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
        
        trig_ampl_i = return_sig_on_trig(t, ampl_i, c_trig_times, baseline_mode)
        
        all_trig_ampl.append(trig_ampl_i.reshape(-1,))
        
    return(np.array(all_trig_ampl))
        
def in_artefact_in_window_all_channels(trig_time,run,rcs):

    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        
        if is_artifact_in_window(anasig_i, trig_time+t_win_start, trig_time+t_win_stop):
            
            print '   SKIP : trig time %f is in artifact window for signal %d '%(trig_time, i)
            return True
            
    return False
    
def clean_trig_all_channels(trig_times, run,rcs):

    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        
        ## cleaning (remove artefact trigs)
        c_trig_times  = clean_trig(trig_times, t_win_start,t_win_stop, anasig_i)
        
        if len(c_trig_times) != len(trig_times):
        
            print "Found artefct for rc %d, removing trigs"%(rcs[i].id)
            trig_times = c_trig_times
        
    return(trig_times)
    

def in_trig_times(time, trig_times):

    #print trig_times
    
    starts = np.array(trig_times) + t_win_start - time
    stops = np.array(trig_times) + t_win_stop - time
    
    #print starts
    #print stops
        
    ok1  = np.any(((starts<=t_win_start)&(stops>t_win_start)) | ((starts<=t_win_stop)&(stops>t_win_stop)))
    # or contain
    ok2 = np.any((starts>=t_win_start ) & (stops<t_win_stop ))
    
    return ok1 or ok2
    
    
def return_random_trig_ampl(run,rcs,trial_type = 'all'):
    
    if trial_type == 'all':
        trig_times = [ trial.triggered_odor_time for trial in run.trials ]
        
    else:
        trig_times = [ trial.triggered_odor_time for trial in run.trials if trial.score_recognition == trial_type ]
    
    c_trig_times = clean_trig_all_channels(trig_times,run,rcs)
    
    print len(c_trig_times)
    
    fake_trig_times = []
    
    count = 0
    
    anasig_0 = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[0].id)[0]
    t0, ampl_0 = get_envelope_in_band(anasig_0.id, f_start, f_stop, method = envelop_method)
    
    while len(fake_trig_times) < len(c_trig_times) and count < 1000:
        
        fake_trig_time =  np.random.rand()*(t0[-1]-t0[0]-t_win_stop+t_win_start)-t_win_start+t0[0]
    
        print fake_trig_time
        
        if not in_artefact_in_window_all_channels(fake_trig_time,run,rcs) and not in_trig_times(fake_trig_time,c_trig_times):
            
            fake_trig_times.append(fake_trig_time)
            print "Not in artefact window and out of trig times"
            
        count += 1
                
    if count == 1000:
        
        print "Warning, after 1000 iterations, could not found time neither in trig epochs nor in artefact windows"
        
        exit()
        
        
    rand_trig_ampl = []
    
    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        t, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
        
        trig_ampl_i = return_sig_on_trig(t, ampl_i, fake_trig_times, baseline_mode)
        
        rand_trig_ampl.append(trig_ampl_i.reshape(-1,))
        
    return(np.array(rand_trig_ampl))
    
def return_all_ampl_by_trig(run,rcs,trial_type = 'all'):
    
    if trial_type == 'all':
        trig_times = [ trial.triggered_odor_time for trial in run.trials ]
        
    else:
        trig_times = [ trial.triggered_odor_time for trial in run.trials if trial.score_recognition == trial_type ]
    
    print trig_times
    print len(trig_times)
    
    c_trig_times = clean_trig_all_channels(trig_times,run,rcs)
    
    print c_trig_times
    print len(c_trig_times)
    
    all_trig_ampl = []
    
    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        t, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
        
        ## cleaning (remove artefact trigs)
        trig_ampl_i = return_sig_on_trig(t, ampl_i, c_trig_times, baseline_mode)
        
         ## no cleaning (all trigs)
        #trig_ampl_i = return_sig_on_trig(t, ampl_i, trig_times, baseline_mode)
        
        all_trig_ampl.append(trig_ampl_i)
        
    return(np.array(all_trig_ampl))
        
        
        
def return_random_ampl_by_trig(run,rcs):
    
    anasig_0 = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[0].id)[0]
    t0, ampl_0 = get_envelope_in_band(anasig_0.id, f_start, f_stop, method = envelop_method)
    
    trig_times = [ trial.triggered_odor_time for trial in run.trials ]
    
    c_trig_times = clean_trig_all_channels(trig_times,run,rcs)
    
    print len(c_trig_times)
    
    fake_trig_times = []
    
    count = 0
    
    while len(fake_trig_times) < len(c_trig_times) and count < 1000:
        
        fake_trig_time =  np.random.rand()*(t0[-1]-t0[0]-t_win_stop+t_win_start)-t_win_start+t0[0]
    
        print fake_trig_time
        
        if not in_artefact_in_window_all_channels(fake_trig_time,run,rcs) and not in_trig_times(fake_trig_time,c_trig_times):
            
            fake_trig_times.append(fake_trig_time)
            print "Not in artefact window and out of trig times"
            
        count += 1
                
    if count == 1000:
        
        print "Warning, after 1000 iterations, could not found time neither in trig epochs nor in artefact windows"
        
        exit()
        
        
    rand_trig_ampl = []
    
    for i in range(len(rcs)):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        t, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
        
        trig_ampl_i = return_sig_on_trig(t, ampl_i, fake_trig_times, baseline_mode)
        
        rand_trig_ampl.append(trig_ampl_i)
        
    return(np.array(rand_trig_ampl))
    
    
    
def compute_cormat(np_all_ampl):

    print np_all_ampl.shape
    
    n = np_all_ampl.shape[0]

    cormat = np.zeros((n,n))

    for (i, j) in itertools.combinations(range(n), 2):
        
        r, p_val = scipy.stats.pearsonr(np_all_ampl[i,], np_all_ampl[j,])
        
        cormat [j,i] = cormat[i,j] = r
        
        print i,j,r
        
    print cormat.shape
    
    return cormat
##########################################################################################################################################################################################################################

def inband_amplitude_correlation_full_sig():
    
    subject = session.query(Subject).filter_by(name = 'FERJ').first()
    run = subject.runs.filter_by(exp= 'E').filter_by(index=1).first()
    #~ print run.subject
    #~ print run
    rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    rcs = [ ]
    for rcg in rcgs:
        rcs.extend(rcg.recordingchannels)
    
    #rcs = rcs[:40]
    rc_names = [rc.name for rc in rcs ]
    n = len(rcs)
    
    ################# full sig (removing artefact and trig windows)
    ts_correct_non_trig_sigs_file = os.path.join(correl_analysis_path,'ts_correct_non_trig_sigs.npy')
    ts_correct_non_trig_ampl_file = os.path.join(correl_analysis_path,'ts_correct_non_trig_ampl.npy')
        
    if not os.path.exists(ts_correct_non_trig_ampl_file) or not os.path.exists(ts_correct_non_trig_sigs_file):
            
        np_correct_non_trig_sigs,np_correct_non_trig_ampl = return_correct_non_trig_ampl_sigs(run,rcs)
        
        np.save(ts_correct_non_trig_ampl_file,np_correct_non_trig_ampl)
        np.save(ts_correct_non_trig_sigs_file,np_correct_non_trig_sigs)
        
    else:
        np_correct_non_trig_ampl = np.load(ts_correct_non_trig_ampl_file)
        np_correct_non_trig_sigs = np.load(ts_correct_non_trig_sigs_file)
        
    print np_correct_non_trig_ampl.shape
    
    #### correl correct ampl 
    cormat_correct_non_trig_ampl_file =  os.path.join(correl_analysis_path,'cormat_correct_non_trig_ampl.npy')
    
    if not os.path.exists(cormat_correct_non_trig_ampl_file):
        
        cormat_correct_non_trig_ampl = compute_cormat(np_correct_non_trig_ampl)
        
        np.save(cormat_correct_non_trig_ampl_file,cormat_correct_non_trig_ampl)
        
    else:
        
        print "loading cormat_correct_non_trig_ampl_file"
        
        cormat_correct_non_trig_ampl = np.load(cormat_correct_non_trig_ampl_file)
        
    ## plotting cor mat
    plot_heatmap_cormat_file =  os.path.join(correl_analysis_path,'heatmap_cormat_correct_non_trig_ampl.eps')
    
    plot_ranged_cormat(plot_heatmap_cormat_file,cormat_correct_non_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
   
   #### correl correct sigs
    cormat_correct_non_trig_sigs_file =  os.path.join(correl_analysis_path,'cormat_correct_non_trig_sigs.npy')
    
    if not os.path.exists(cormat_correct_non_trig_sigs_file):
        
        cormat_correct_non_trig_sigs = compute_cormat(np_correct_non_trig_sigs)
        
        np.save(cormat_correct_non_trig_sigs_file,cormat_correct_non_trig_sigs)
        
    else:
        
        print "loading cormat_correct_non_trig_sigs_file"
        
        cormat_correct_non_trig_sigs = np.load(cormat_correct_non_trig_sigs_file)
        
    ## plotting cor mat
    plot_heatmap_cormat_file =  os.path.join(correl_analysis_path,'heatmap_cormat_correct_non_trig_sigs.eps')
    
    plot_ranged_cormat(plot_heatmap_cormat_file,cormat_correct_non_trig_sigs,list_labels = rc_names,fix_full_range = [-1.0,1.0])
   
   
   
def export_inband_trig_amplitude(export_dir,name,exp,index,null_model = 'random_trig'):
    
    subject = session.query(Subject).filter_by(name = name).first()
    run = subject.runs.filter_by(exp= exp).filter_by(index=index).first()
    #~ print run.subject
    #~ print run
    rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    rcs = [ ]
    for rcg in rcgs:
        rcs.extend(rcg.recordingchannels)
    
    #rcs = rcs[:40]
    rc_names = [rc.name for rc in rcs ]
    n = len(rcs)
       
    print n
        
    dir_sess = name + '_' + exp + str(index) 
    
    
    try :
        
        os.makedirs(os.path.join(export_dir,dir_sess))
        
    except OSError:
        print "Warning OS error, dir_sess %s already exists"%(dir_sess)

    ###### np series by trigs
    for trial_type in trial_R_conds:
    #for trial_type in ['cr']:
        
        
        ### concat all trigs in one time series
        print trial_type
        
        #print "trig_ampl"
        #ts_all_trig_ampl_file =  os.path.join(export_dir,dir_sess,'ts_' + trial_type+'_trig_ampl.npy')
        
        #if not os.path.exists(ts_all_trig_ampl_file):
            
            #np_all_trig_ampl = return_all_trig_ampl(run,rcs,trial_type = trial_type)

            #np.save(ts_all_trig_ampl_file,np_all_trig_ampl)
        #else:
            
            #np_all_trig_ampl = np.load(ts_all_trig_ampl_file)
        
        #print np_all_trig_ampl.shape
        
        #### keep by trig
        #print "ampl_by_trig"
        
        #ts_all_ampl_by_trig_file =  os.path.join(export_dir,dir_sess,'ts_' + trial_type+'_ampl_by_trig.npy')
        
        #if not os.path.exists(ts_all_ampl_by_trig_file):
            
            #np_all_ampl_by_trig = return_all_ampl_by_trig(run,rcs,trial_type = trial_type)

            #np.save(ts_all_ampl_by_trig_file,np_all_ampl_by_trig)
        #else:
            
            #np_all_ampl_by_trig = np.load(ts_all_ampl_by_trig_file)
        
        #print np_all_ampl_by_trig.shape
        
        
        ### keep by trig
        print "ampl_by_trial"
        
        ## odor + rest
        return_correct_ampl_by_trig_duration_and_rest(os.path.join(export_dir,dir_sess),run,rcs,trial_type = trial_type)
        
        ## only odor
        #return_correct_ampl_by_trig_duration(os.path.join(export_dir,dir_sess),trial_type = trial_type)
            
        
    channel_names_file = os.path.join(export_dir,dir_sess,'channel_names.txt')
    
    np.savetxt(channel_names_file,rc_names, fmt = '%s')
    
    
    
#def inband_trig_amplitude_correlation(name,exp,index,trial_type,null_model = 'random_trig'):
    
    #subject = session.query(Subject).filter_by(name = name).first()
    #run = subject.runs.filter_by(exp= exp).filter_by(index=index).first()
    ##~ print run.subject
    ##~ print run
    #rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    #rcs = [ ]
    #for rcg in rcgs:
        #rcs.extend(rcg.recordingchannels)
    
    ##rcs = rcs[:40]
    #rc_names = [rc.name for rc in rcs ]
    #n = len(rcs)
        
    #dir_sess = name + '_' + exp + str(index) + '_' + trial_type
    
    
    #try :
        
        #os.makedirs(os.path.join(correl_analysis_path,dir_sess))
        
    #except OSError:
        #print "Warning OS error, dir_sess %s already exists"%(dir_sess)

        
    
    ####### all trigs
    
    #ts_all_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'ts_all_trig_ampl.npy')
    
    #if not os.path.exists(ts_all_trig_ampl_file):
        
        #np_all_trig_ampl = return_all_trig_ampl(run,rcs,trial_type = trial_type)

        #print np_all_trig_ampl.shape
        
        #np.save(ts_all_trig_ampl_file,np_all_trig_ampl)
    #else:
        
        #np_all_trig_ampl = np.load(ts_all_trig_ampl_file)
    
    #print np_all_trig_ampl.shape
    
    
    
    
    
    #cormat_all_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'cormat_all_trig_ampl.npy')
    
    #if not os.path.exists(cormat_all_trig_ampl_file):
        
        #cormat_all_trig_ampl = compute_cormat(np_all_trig_ampl)
        
        #np.save(cormat_all_trig_ampl_file,cormat_all_trig_ampl)
        
    #else:
        
        #print "loading cormat_all_trig_ampl_file"
        
        #cormat_all_trig_ampl = np.load(cormat_all_trig_ampl_file)
        
    ### plotting cor mat
    #plot_heatmap_cormat_all_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'heatmap_cormat_all_trig_ampl.eps')
    
    #plot_ranged_cormat(plot_heatmap_cormat_all_trig_ampl_file,cormat_all_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
   
    ######## random trigs
    
    #if null_model == 'random_trig':
        
        #ts_random_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'ts_random_trig_ampl.npy')
        
        #if not os.path.exists(ts_random_trig_ampl_file):
            
            #np_random_trig_ampl = return_random_trig_ampl(run,rcs,trial_type = trial_type)

            #print np_random_trig_ampl.shape
            
            #np.save(ts_random_trig_ampl_file,np_random_trig_ampl)
        #else:
            
            #np_random_trig_ampl = np.load(ts_random_trig_ampl_file)
        
        #print np_random_trig_ampl.shape
        
        #cormat_random_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'cormat_random_trig_ampl.npy')
        
        #if not os.path.exists(cormat_random_trig_ampl_file):
            
            #cormat_random_trig_ampl = compute_cormat(np_random_trig_ampl)
            
            #np.save(cormat_random_trig_ampl_file,cormat_random_trig_ampl)
            
        #else:
            
            #print "loading cormat_random_trig_ampl_file"
            
            #cormat_random_trig_ampl = np.load(cormat_random_trig_ampl_file)
            
        ### plotting cor mat
        #plot_heatmap_cormat_random_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'heatmap_cormat_random_trig_ampl.eps')
        
        #plot_ranged_cormat(plot_heatmap_cormat_random_trig_ampl_file,cormat_random_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
    
    
    
        ### plotting diff cor mat
        #cormat_diff_trig_ampl = cormat_all_trig_ampl - cormat_random_trig_ampl
        
        #plot_heatmap_cormat_diff_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'heatmap_cormat_diff_trig_ampl.eps')
        
        #plot_ranged_cormat(plot_heatmap_cormat_diff_trig_ampl_file,cormat_diff_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
    
    ##elif null_model == 'correct_non_trig_ampl':
            
        ##ts_all_correct_ampl_file = os.path.join(correl_analysis_path,dir_sess,'ts_correct_non_trig_ampl.npy')
            
        ##if not os.path.exists(ts_all_correct_ampl_file):
                
            ##_,np_all_ampl = return_correct_non_trig_ampl_sigs(run,rcs)
            
            ##np.save(ts_all_correct_ampl_file,np_all_ampl)
            
        ##else:
            ##np_all_ampl = np.load(ts_all_correct_ampl_file)
            
        ##print np_all_ampl.shape
        
        ###### correl correct ampl 
        ##cormat_correct_non_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'cormat_correct_non_trig_ampl.npy')
        
        ##if not os.path.exists(cormat_correct_non_trig_ampl_file):
            
            ##cormat_correct_non_trig = compute_cormat(np_all_ampl)
            
            ##np.save(cormat_correct_non_trig_ampl_file,cormat_correct_non_trig_ampl)
            
        ##else:
            
            ##print "loading cormat_correct_non_trig_ampl_file"
            
            ##cormat_correct_non_trig_ampl = np.load(cormat_correct_non_trig_ampl_file)
            
        #### plotting cor mat
        ##plot_heatmap_cormat_file =  os.path.join(correl_analysis_path,dir_sess,'heatmap_cormat_correct_non_trig_ampl.eps')
        
        ##plot_ranged_cormat(plot_heatmap_cormat_file,cormat_correct_non_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
    
    
    
        #### plotting diff cor mat
        ##cormat_diff_trig_ampl_correct_non_trig_ampl = cormat_all_trig_ampl - cormat_correct_non_trig_ampl
        
        ##plot_heatmap_cormat_diff_trig_ampl_correct_non_trig_ampl_file =  os.path.join(correl_analysis_path,dir_sess,'heatmap_cormat_diff_trig_ampl_correct_non_trig_ampl.eps')
        
        ##plot_ranged_cormat(plot_heatmap_cormat_diff_trig_ampl_correct_non_trig_ampl_file,cormat_diff_trig_ampl_correct_non_trig_ampl,list_labels = rc_names,fix_full_range = [-1.0,1.0])
    
   
#def inband_amplitude_correlation_by_trig(name,exp,index,trial_type):
    
    #subject = session.query(Subject).filter_by(name = name).first()
    #run = subject.runs.filter_by(exp= exp).filter_by(index=index).first()
    ##~ print run.subject
    ##~ print run
    #rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    #rcs = [ ]
    #for rcg in rcgs:
        #rcs.extend(rcg.recordingchannels)
    
    ##rcs = rcs[:40]
    #rc_names = [rc.name for rc in rcs ]
    #n = len(rcs)
        
    #np_all_ampl_by_trig = return_all_ampl_by_trig(run,rcs,trial_type = trial_type)

    #print np_all_ampl_by_trig
    #print np_all_ampl_by_trig.shape
    
    #np_rand_ampl_by_trig = return_random_ampl_by_trig(run,rcs)

    #print np_rand_ampl_by_trig
    #print np_rand_ampl_by_trig.shape
    
    #list_diff = []
    
    
    #for (i, j) in itertools.combinations(range(n), 2):
        
        #Z_values = []
        #Z_rand_values = []
        
        #for trial in range(np_all_ampl_by_trig.shape[1]):
        
            #r, p_val = scipy.stats.pearsonr(np_all_ampl_by_trig[i,trial,], np_all_ampl_by_trig[j,trial,])
        
            ##print r,p_val
            
            #Z_values.append(np.arctanh(r))
            
            #r_rand, p_val_rand = scipy.stats.pearsonr(np_rand_ampl_by_trig[i,trial,], np_rand_ampl_by_trig[j,trial,])
        
            ##print r_rand, p_val_rand
            
            #Z_rand_values.append(np.arctanh(r_rand))
            
        ##print Z_values
        
        
        
        #### test if both samples are different
        #t,p_val = ttest_ind(Z_values,Z_rand_values)
        
        
        #list_diff.append([i,j,p_val,np.sign(t)])
        
    #np_list_diff = np.array(list_diff)
   
    #signif_code = return_signif_code(np_list_diff[:,2],uncor_alpha = 0.001,fdr_alpha = 0.05, bon_alpha = 0.05)
        
    #print signif_code
        
    #print np.sum(signif_code == 0.0),np.sum(signif_code == 1.0),np.sum(signif_code == 2.0),np.sum(signif_code == 3.0),np.sum(signif_code == 4.0)
    
    #np_list_diff[:,3] = np_list_diff[:,3] * signif_code
    
    #print np.sum(np_list_diff[:,3] == 0.0)
    #print np.sum(np_list_diff[:,3] == 1.0),np.sum(np_list_diff[:,3] == 2.0),np.sum(np_list_diff[:,3] == 3.0),np.sum(np_list_diff[:,3] == 4.0)
    #print np.sum(np_list_diff[:,3] == -1.0),np.sum(np_list_diff[:,3] == -2.0),np.sum(np_list_diff[:,3] == -3.0),np.sum(np_list_diff[:,3] == -4.0)
    
    
    
    
    #signif_signed_adj_mat = np.zeros((n,n),dtype = 'int')
    
    #signif_i = np.array(np_list_diff[:,0],dtype = int)
    #signif_j = np.array(np_list_diff[:,1],dtype = int)
    
    #signif_sign = np.array(np_list_diff[:,3],dtype = int)
    
    #print signif_i,signif_j
    
    #print signif_signed_adj_mat[signif_i,signif_j] 
    
    ##print signif_sign
    
    
    
    #signif_signed_adj_mat[signif_i,signif_j] = signif_signed_adj_mat[signif_j,signif_i] = signif_sign
    
    #print signif_signed_adj_mat
            
            
    ### plotting diff cor mat
    
    #plot_heatmap_signif_signed_adj_mat_by_trig_file =  os.path.join(correl_analysis_path,'heatmap_signif_signed_adj_mat_by_trig.eps')
    
    #plot_int_mat(plot_heatmap_signif_signed_adj_mat_by_trig_file,signif_signed_adj_mat,list_labels = rc_names,fix_full_range = [-4.0,4.0])
   
def export_ts_retrieval_all_cond():

    for subj_id in subject_ids:
    #for subj_id in ['SEMC']:
        
        for R_index in ['1','2','3']:
        
            export_inband_trig_amplitude(export_dir = base_dir,name = subj_id,exp = 'R',index = R_index)
            
    #export_inband_trig_amplitude(export_dir = base_dir,name = 'FERJ',exp = 'R',index = '1')

if __name__ == '__main__':
    #inband_amplitude_correlation_full_sig()
    #inband_trig_amplitude_correlation(name = 'FERJ',exp = 'R',index= '1',trial_type = 'hit',null_model = 'random_trig')
    #inband_amplitude_correlation_by_trig(name = 'FERJ',exp = 'E',index= '1',trial_type = 'all')
    
    
    ### mise en place pipeline:
    export_ts_retrieval_all_cond()