# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 

from connection import *

from params import joblib_cachedir
from params import subject_ids,trial_R_conds,trial_E_conds

from params import envelop_method,baseline_mode

from params import t_win_start_odor,t_win_stop_odor
from params import t_win_start_rest,t_win_stop_rest

from OpenElectrophy.timefrequency import TimeFreq

#from inband_amplitude_around_trigger import get_envelope_in_band

from artifact_detection import is_artifact_in_window,detect_artifact_on_one_sig,return_artifact_timings,clean_trig

from timefreqtools.computations import get_timefreq,get_envelope_in_band

from timefreqtools.computations import get_bipolar_timefreq,get_bipolar_envelope_in_band

import timefreqtools
timefreqtools.set_session(session, dbinfo)

#from dmgraphanalysis_nodes.utils_stats import return_signif_code

#from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat,plot_int_mat
        

        
  
import itertools

#~ from artifact_detection import is_artifact_in_window

import scipy.signal
from scipy.stats import ttest_ind
from matplotlib import pyplot



#### compute envelope amplitude in frequency band with hilbert or tfr 
def compute_ampl_sigs(run,rcs,f_start, f_stop):
    
    ampl_sigs = []
    
    if len(rcs) != 0:
        
        for i in range(len(rcs)):
            
            print i
            
            anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
            
            times, ampl_i = get_envelope_in_band(anasig_i.id, f_start, f_stop, method = envelop_method)
            
            ampl_sigs.append(ampl_i)
            
    else:
        
        times = np.array([])
    #print sigs
    
    return np.array(ampl_sigs,dtype = float),times
    
def compute_neo_sigs(run,rcs):
    
    neo_sigs = []
    
    if len(rcs) != 0:
        
        for i in range(len(rcs)):
            
            print i
            
            anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
            
            anasig = session.query(AnalogSignal).get(anasig_i.id)
            
            neosig = anasig.to_neo()
            
            times = neosig.times.rescale('s').magnitude
            
            neo_sigs.append(neosig)
            
    else:
        
        times = np.array([])
    #print sigs
    
    return neo_sigs,times
    
def compute_bip_ampl_sigs(run,bip_rc_ids,f_start, f_stop):
    
    bip_ampl_sigs = []
    
    if len(bip_rc_ids) != 0:
        
        #for chan_name in [bip_rc_ids[0]]:
        for chan_name in bip_rc_ids:
            
            print chan_name
            
            rc_i,rc_j = chan_name.split("-")
            
            print rc_i,rc_j
            
            anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc_i)[0]
            
            anasig_j = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc_j)[0]
                        
            print anasig_i.id,anasig_j.id
            
            times, ampl_i = get_bipolar_envelope_in_band(anasig_i.id,anasig_j.id, f_start, f_stop)
            
            bip_ampl_sigs.append(ampl_i)
            
    else:
        
        times = np.array([])
    #print sigs
    
    return np.array(bip_ampl_sigs,dtype = float),times
    
    
def return_bip_ampl_sigs(cur_dir,run,bip_rc_ids,f_start, f_stop):
    
    np_bip_ampl_sigs_file = os.path.join(cur_dir,'np_bip_ampl_sigs.npy')
    
    timings_file = os.path.join(cur_dir,'timings.npy')
        
    if not os.path.exists(np_bip_ampl_sigs_file) or not os.path.exists(timings_file):
        
        np_bip_ampl_sigs,timings = compute_bip_ampl_sigs(run,bip_rc_ids,f_start, f_stop)
        
        np.save(np_bip_ampl_sigs_file,np_bip_ampl_sigs)
        
        np.save(timings_file,timings)
    
    else:
        
        np_bip_ampl_sigs = np.load(np_bip_ampl_sigs_file)
        
        timings = np.load(timings_file)
        
    return np_bip_ampl_sigs,timings
    
def return_ampl_sigs(cur_dir,run,rcs,f_start, f_stop):
    
    np_ampl_sigs_file = os.path.join(cur_dir,'np_ampl_sigs.npy')
    
    timings_file = os.path.join(cur_dir,'timings.npy')
        
    if not os.path.exists(np_ampl_sigs_file) or not os.path.exists(timings_file):
        
        np_ampl_sigs,timings = compute_ampl_sigs(run,rcs,f_start, f_stop)
        
        np.save(np_ampl_sigs_file,np_ampl_sigs)
        
        np.save(timings_file,timings)
    
    else:
        
        np_ampl_sigs = np.load(np_ampl_sigs_file)
        
        timings = np.load(timings_file)
        
    return np_ampl_sigs,timings
    
############################ timings ########################################
def compute_timings(run,rcs):

    times = np.array([])
    
    for i in range(len(rcs)):
        
        #print i
        
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
        
        neosig_i = anasig_i.to_neo()
        
        if times.shape[0] < neosig_i.times.rescale('s').magnitude.shape[0]:
        
            times = neosig_i.times.rescale('s').magnitude
        
        #print times
        
    return times
        

def return_timings(timings_file,run,rcs):

    if not os.path.exists(timings_file):
        
        timings = compute_timings(run,rcs)
        
        np.save(timings_file,timings)
                        
        
    else:
        
        timings = np.load(timings_file)
        
    return timings

######################### electrodes names and coords     

def return_bip_electrode_names_and_coords(rcs,cur_dir):
    
    ########## channel names ###########
    
    rc_names = [rc.name for rc in rcs ]
    
    bip_rc_names = ["-".join(a) for a in zip(rc_names[:-1],rc_names[1:])]
    
    print len(bip_rc_names)
    
    channel_names_file = os.path.join(cur_dir,'bip_channel_names.txt')
    
    np.savetxt(channel_names_file,bip_rc_names, fmt = '%s')
        
    ########### channel ids ############
    rc_ids = [rc.id for rc in rcs ]
    
    print rc_ids
    
    bip_rc_ids = ["-".join(a) for a in zip(map(str,rc_ids[:-1]),map(str,rc_ids[1:]))]
    
    print len(bip_rc_ids)
    
    channel_ids_file = os.path.join(cur_dir,'bip_channel_ids.txt')
    
    np.savetxt(channel_ids_file,bip_rc_ids, fmt = '%s')
    
    ########## channel coords ##########
    
    channel_coords_file = os.path.join(cur_dir,'bip_channel_coords.txt')
    
    rc_coords = []
    for rc in rcs:
        
        #print rc.coordinate
        
        if rc.coordinate is None:
        
            print "warning, coordiantes is type None"
            a = np.empty(shape = (3), dtype = float)
            a.fill(np.nan)
            
            rc_coords.append(a)
            
            #print a
            #print rc_coords
            
        else:
            rc_coords.append(rc.coordinate.magnitude)
            
    #print rc_coords
    
    np_rc_coords = np.array(rc_coords)
    
    
    bip_np_rc_coords = (np_rc_coords[:-1,] + np_rc_coords[1:,]) / 2
    
    print bip_np_rc_coords.shape
    
    bip_rc_coords = bip_np_rc_coords.tolist()
    
    #### if everything works well... (all coordinates are filled)
    #rc_coords = [rc.coordinate.magnitude for rc in rcs if rc.coordinate is not 'None']
    
    np.savetxt(channel_coords_file,bip_rc_coords, fmt = '%.2f')
    
    print len(bip_rc_names),len(bip_rc_coords)
    
    if len(bip_rc_names) != len(bip_rc_coords):
    
        print "$$$$$$$$$$$$$$$$$$ Warning, descrepancy between names and coords"
        
        sys.exit()
        
    return bip_rc_names,bip_rc_coords,bip_rc_ids
    
    
def return_electrode_names_and_coords(rcs,cur_dir):
    
    ########## channel names ###########
    rc_names = [rc.name for rc in rcs ]
    
    rc_descriptions = [rc.description for rc in rcs]
            
    channel_names_file = os.path.join(cur_dir,'channel_names.txt')
    
    np.savetxt(channel_names_file,rc_names, fmt = '%s')
    
    ########## channel coords ##########
    
    channel_coords_file = os.path.join(cur_dir,'channel_coords.txt')
    
    rc_coords = []
    
    ######### channel descriptions ############
    
    channel_descriptions_file = os.path.join(cur_dir,'channel_descriptions.txt')
    
    np.savetxt(channel_descriptions_file,rc_descriptions, fmt = '%s')
    
    
    for rc in rcs:
        
        #print rc.coordinate
        
        if rc.coordinate is None:
        
            print "warning, coordiantes is type None"
            a = np.empty(shape = (3), dtype = float)
            a.fill(np.nan)
            
            rc_coords.append(a)
            
            #print a
            #print rc_coords
            
        else:
            rc_coords.append(rc.coordinate.magnitude)
            
    #print rc_coords
    
    #### if everything works well... (all coordinates are filled)
    #rc_coords = [rc.coordinate.magnitude for rc in rcs if rc.coordinate is not 'None']
    
    np.savetxt(channel_coords_file,rc_coords, fmt = '%f')
    
    print len(rc_names),len(rc_coords)
    
    if len(rc_names) != len(rc_coords):
    
        print "$$$$$$$$$$$$$$$$$$ Warning, descrepancy between names and coords"
        
        sys.exit()
        
    return rc_names,rc_coords,rc_descriptions
   
 ############### keep_all_non_art ###########################
 
def compute_keep_from_timings(timings,art_starts_stops):

    keep = np.ones(shape = timings.shape, dtype = 'bool')
    
    for art_seq in art_starts_stops:
        
        keep[(timings >= art_seq[0]) & (timings < art_seq[1])] = False
        
    return keep
    
def compute_non_art(run,rcs, timings):
    
    keep = np.ones(shape = (len(rcs),timings.shape[0]), dtype = 'bool')
    
    #print keep.shape
    
    #### computing artefacts
    for i,rc in enumerate(rcs):
        anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc.id)[0]
        
        art_start,art_stop = return_artifact_timings(anasig_i)
        
        keep[i,:] = compute_keep_from_timings(timings,zip(art_start,art_stop))
        
        #print np.sum(keep[i,:] == 0)
        
    #print np.sum(keep == 0,axis = 1)
    
    return keep, np.array( np.sum(keep,axis = 0) == keep.shape[0],dtype = 'bool')
    
    
def compute_bip_non_art(run,bip_rc_ids, timings):
    
    keep = np.ones(shape = (len(bip_rc_ids),timings.shape[0]), dtype = 'bool')
    
    #print keep.shape
    
    #### computing artefacts
    for i,bip_rc_id in enumerate(bip_rc_ids):
    #for i,bip_rc_id in enumerate([bip_rc_ids[0]]):
        
        print bip_rc_id
        
        rc_id_i,rc_id_j = bip_rc_id.split("-")
        
        print  rc_id_i,rc_id_j
        
        anasig1 = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc_id_i)[0]
        
        art_start1,art_stop1 = return_artifact_timings(anasig1)
        
        keep1 = compute_keep_from_timings(timings,zip(art_start1,art_stop1))
        
        print keep1
        
        print np.sum(keep1 == False)
        
        #print np.where(keep1 == False)
        
        anasig2 = run.segments[0].analogsignals.filter_by(recordingchannel_id = rc_id_j)[0]
        
        art_start2,art_stop2 = return_artifact_timings(anasig2)
        
        keep2 = compute_keep_from_timings(timings,zip(art_start2,art_stop2))
        
        print np.sum(keep2 == False)
        #print np.where(keep2 == False)
        
        print np.where(np.logical_or(keep1 == False,keep2 == False) == True)
        
        #print np.logical_not(keep1 == False,keep2 == False) == True)
        
        keep[i,:] = np.logical_not(np.logical_or(keep1 == False,keep2 == False) == True)
        
        print np.sum(keep[i,:] == False)
        
    #print np.sum(keep == 0,axis = 1)
    
    return keep, np.array( np.sum(keep,axis = 0) == keep.shape[0],dtype = 'bool')
    
    
    
def return_non_art(export_dir,run,rcs,timings):
            
        keep_all_non_art_file = os.path.join(export_dir,'keep_all_non_art.npy')
        
        non_art_file = os.path.join(export_dir,'non_art.npy')
        
        if not os.path.exists(keep_all_non_art_file) or not os.path.exists(non_art_file):
            
            non_art,keep_all_non_art = compute_non_art(run,rcs,timings)
            
            np.save(keep_all_non_art_file,keep_all_non_art)
            
            np.save(non_art_file,non_art)
            
        else:
            
            keep_all_non_art = np.load(keep_all_non_art_file)
            
            non_art = np.load(non_art_file)
            
            
        return non_art,keep_all_non_art
        
def return_bip_non_art(export_dir,run,bip_rc_ids,timings):
            
        keep_all_bip_non_art_file = os.path.join(export_dir,'keep_all_bip_non_art.npy')
        bip_non_art_file = os.path.join(export_dir,'bip_non_art.npy')
        
        if not os.path.exists(keep_all_bip_non_art_file) or not os.path.exists(bip_non_art_file):
            
            bip_non_art,keep_all_bip_non_art = compute_bip_non_art(run,bip_rc_ids,timings)
            
            np.save(keep_all_bip_non_art_file,keep_all_bip_non_art)
            
            np.save(bip_non_art_file,bip_non_art)
            
        else:
            
            keep_all_bip_non_art = np.load(keep_all_bip_non_art_file)
            
            bip_non_art = np.load(bip_non_art_file)
            
            
        return bip_non_art,keep_all_bip_non_art

    ################################################ trigs and rest ####################################################################

def compute_correct_trigs_and_rest(run,timings,keep_all_non_art,trial_type = 'all'):
    
    if trial_type == 'all':
        
        start_odor = [ trial.first_inspiration_time + t_win_start_odor for trial in run.trials ]
        stop_odor = [ trial.first_inspiration_time + t_win_stop_odor for trial in run.trials]
        start_rest = [ trial.first_inspiration_time + t_win_start_rest for trial in run.trials]
        stop_rest = [ trial.first_inspiration_time + t_win_stop_rest for trial in run.trials]
        
    else:
        
        start_odor = [ trial.first_inspiration_time + t_win_start_odor for trial in run.trials if trial.score_recognition == trial_type ]
        stop_odor = [ trial.first_inspiration_time + t_win_stop_odor for trial in run.trials if trial.score_recognition == trial_type ]
        start_rest = [ trial.time + t_win_start_rest for trial in run.trials if trial.score_recognition == trial_type ]
        stop_rest = [ trial.time + t_win_stop_rest for trial in run.trials if trial.score_recognition == trial_type ]
        
    print len(start_odor), len(stop_odor), len(start_rest ), len(stop_rest)
    
    correct_trigs = []
    
    correct_rest = []
    
    nb_removed_trigs_odor = 0

    nb_removed_trigs_rest = 0
    
    for i,odor_range in enumerate(zip(start_odor,stop_odor)):
        
        sel_odor = (timings > odor_range[0]) & (timings <= odor_range[1])
        
        sel_rest = (timings > start_rest[i]) & (timings <= stop_rest[i])
        
        #print keep_all_non_art[sel_odor]
        
        if np.sum(keep_all_non_art[sel_odor]) == np.sum(sel_odor == True):
            
            correct_trigs.append(odor_range)
            
        else:
            
            nb_removed_trigs_odor = nb_removed_trigs_odor + 1

        if np.sum(keep_all_non_art[sel_rest]) == np.sum(sel_rest == True):
        #if np.sum(sel_odor == True)*0.98 <= np.sum(keep_all_non_art[sel_odor]) and np.sum(sel_rest == True)*0.9  <= np.sum(keep_all_non_art[sel_rest]):
        
            correct_rest.append((start_rest[i],stop_rest[i]))
            
        else:
            
            nb_removed_trigs_rest = nb_removed_trigs_rest + 1
        
    print len(correct_trigs),len(correct_rest),nb_removed_trigs_odor, nb_removed_trigs_rest, 
    
    return correct_trigs,correct_rest,nb_removed_trigs_odor, nb_removed_trigs_rest,

def return_sess_correct_trigs_and_rest(run,sess_timings,sess_keep_all_non_art,trial_types):
        
        #### compute all correct trigs and rest periods 
        sess_correct_trigs = []
        
        sess_correct_rest = []

        nb_removed_trigs_odor_by_trial = []

        nb_removed_trigs_rest_by_trial = []
        
        for trial_type in trial_types:
        #for trial_type in ['hit']:
            
            print trial_type
            
            ## odor + rest
            correct_trigs,correct_rest,nb_removed_trigs_odor, nb_removed_trigs_rest = compute_correct_trigs_and_rest(run,sess_timings,sess_keep_all_non_art,trial_type = trial_type)
            
            sess_correct_trigs.append(correct_trigs)
            
            sess_correct_rest.append(correct_rest)
        
            nb_removed_trigs_odor_by_trial.append(nb_removed_trigs_odor)
            nb_removed_trigs_rest_by_trial.append(nb_removed_trigs_rest)
            
        print sess_correct_trigs
        
        return sess_correct_trigs,sess_correct_rest,nb_removed_trigs_odor_by_trial, nb_removed_trigs_rest_by_trial
        
def compute_trial_ampl_sigs(correct_ampl_sigs,correct_trigs,sess_timings):

    trial_ampl_sigs = []
    
    print correct_trigs
    
    for trig in correct_trigs:
        
        print trig
        
        sel = (trig[0] < sess_timings) & (sess_timings <= trig[1])
        
        trial_ampl_sigs.append(correct_ampl_sigs[:,sel])
        
    return np.array(trial_ampl_sigs)

    
    
def return_sess_correct_ampl_trigs(write_path,sess_correct_trigs,sess_correct_ampl_sigs,sess_timings,trial_types):

    print sess_correct_ampl_sigs.shape
    
    sess_correct_ampl_trigs = []
    
    for i,trial_type in enumerate(trial_types):
        
        odor_ampl_file = os.path.join(write_path,'correct_ts_' + trial_type+'_ampl_by_odor_trigs.npy')
        
        if not os.path.exists(odor_ampl_file):
            
            trial_ampl_sigs = compute_trial_ampl_sigs(sess_correct_ampl_sigs,sess_correct_trigs[i],sess_timings)
            
            print trial_ampl_sigs.shape
            
            np.save(odor_ampl_file,trial_ampl_sigs)
        else:
            
            trial_ampl_sigs = np.load(odor_ampl_file)
            
        sess_correct_ampl_trigs.append(trial_ampl_sigs)
        
    return sess_correct_ampl_trigs
            
        
def return_sess_correct_ampl_rest(write_path,sess_correct_trigs,sess_correct_ampl_sigs,sess_timings,trial_types):

    print sess_correct_ampl_sigs.shape
    
    sess_correct_ampl_trigs = []
    
    for i,trial_type in enumerate(trial_types):
        
        rest_ampl_file = os.path.join(write_path,'correct_ts_' + trial_type+'_ampl_by_rest_trigs.npy')
        
        if not os.path.exists(rest_ampl_file):
            
            trial_ampl_sigs = compute_trial_ampl_sigs(sess_correct_ampl_sigs,sess_correct_trigs[i],sess_timings)
            
            print trial_ampl_sigs.shape
            
            np.save(rest_ampl_file,trial_ampl_sigs)
        else:
            
            trial_ampl_sigs = np.load(rest_ampl_file)
            
        sess_correct_ampl_trigs.append(trial_ampl_sigs)
        
    return sess_correct_ampl_trigs
    
############################################################################## proba density #############################################################################################

def compute_trial_proba_density(correct_neo_sigs,correct_trigs,f_start,f_stop,concatenate_chunks = True):

    from timefreqtools.computations import get_neo_timefreq_proba_density

    print len(correct_neo_sigs)
    print correct_trigs
    
    all_proba_density = []
    
    if concatenate_chunks == True:
            
        for i,correct_neo_sig in enumerate(correct_neo_sigs):
        
            print i
            
            dist, freqs, bins, m, iso_curves = get_neo_timefreq_proba_density(neosig = correct_neo_sig,f_start=f_start, f_stop=f_stop, nb_line = int(f_stop - f_start),ylim = 50.,chunks = correct_trigs)
            
            print dist.shape
            
            all_proba_density.append(dist)
            
        return np.array(all_proba_density,dtype = float)
        
    else: 
    
        for i,correct_neo_sig in enumerate(correct_neo_sigs):
        
            print i
            
            list_dist = get_neo_timefreq_proba_density(neosig = correct_neo_sig,f_start=f_start, f_stop=f_stop, nb_line = int(f_stop - f_start),ylim = 50.,chunks = correct_trigs,concatenate_chunks = False)
            
            print len(list_dist)
            
            all_proba_density.append(list_dist)
            
        return np.array(all_proba_density,dtype = float)
        
def return_sess_proba_density_odor(write_path,sess_correct_trigs,sess_correct_neo_sigs,trial_types,f_start,f_stop):

    for i,trial_type in enumerate(trial_types):
        
        odor_proba_density_file = os.path.join(write_path,'proba_density_' + trial_type + '_by_odor_trigs.npy')
        
        if not os.path.exists(odor_proba_density_file):
            
            proba_density_trigs = compute_trial_proba_density(sess_correct_neo_sigs,sess_correct_trigs[i],f_start,f_stop)
            
            print proba_density_trigs.shape
            
            np.save(odor_proba_density_file,proba_density_trigs)
            
            
        odor_chunks_proba_density_file = os.path.join(write_path,'proba_density_chunks_' + trial_type + '_by_odor_trigs.npy')
        
        if not os.path.exists(odor_chunks_proba_density_file):
            
            proba_density_chunks_trigs = compute_trial_proba_density(sess_correct_neo_sigs,sess_correct_trigs[i],f_start,f_stop,concatenate_chunks = False)
            
            print proba_density_chunks_trigs.shape
            
            np.save(odor_chunks_proba_density_file,proba_density_chunks_trigs)
            
def return_sess_proba_density_rest(write_path,sess_correct_trigs,sess_correct_neo_sigs,trial_types,f_start,f_stop):

    for i,trial_type in enumerate(trial_types):
        
        ### concat all signals
        rest_proba_density_file = os.path.join(write_path,'proba_density_' + trial_type+'_by_rest_trigs.npy')
        
        if not os.path.exists(rest_proba_density_file):
            
            proba_density_trigs = compute_trial_proba_density(sess_correct_neo_sigs,sess_correct_trigs[i],f_start,f_stop)
            
            print proba_density_trigs.shape
            
            np.save(rest_proba_density_file,proba_density_trigs)
        
        ### by_chunks
        rest_chunks_proba_density_file = os.path.join(write_path,'proba_density_chunks_' + trial_type + '_by_rest_trigs.npy')
        
        if not os.path.exists(rest_chunks_proba_density_file):
            
            proba_density_chunks_trigs = compute_trial_proba_density(sess_correct_neo_sigs,sess_correct_trigs[i],f_start,f_stop,concatenate_chunks = False)
            
            print proba_density_chunks_trigs.shape
            
            np.save(rest_chunks_proba_density_file,proba_density_chunks_trigs)
            
        
##########################################################################################################################################################################################################################
def compute_subj_correct_bip_ampl(export_dir, f_start, f_stop, name = 'SEMC',exp = 'R',thr_nb_artefact = 120):

    subject = session.query(Subject).filter_by(name = name).first()
    
    rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    rcs = [ ]
    for rcg in rcgs:
        rcs.extend(rcg.recordingchannels)
    
    
    
    dir_subj = os.path.join(export_dir,name)
    
    log_file = os.path.join(dir_subj,'log.txt')
    
    
    
    try :
        print "Creating dir_subj %s"%(dir_subj)
        
        os.makedirs(dir_subj)
        
    except OSError:
        print "Warning OS error, dir_subj %s already exists"%(dir_subj)


    #### electrode names and electrodes
    bip_rc_names,bip_rc_coords,bip_rc_ids = return_bip_electrode_names_and_coords(rcs,dir_subj)

    print len(bip_rc_names),len(bip_rc_coords),len(bip_rc_ids)
    
    
    #### running
    list_bip_ampl_sigs = []
    
    list_sess_timings = []
    
    list_keep_all_bip_non_art = []
    
    list_bip_non_art = []
    
    if exp == 'R':
    
        sess_indexes = ['1','2','3']
        
        trial_types = trial_R_conds
        
    elif exp == 'E':
        
        sess_indexes = ['2']
        
        trial_types = trial_E_conds
        
    
    for sess_index in sess_indexes:
        
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
            
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
        
        try :
            print "Creating dir_sess %s"%(dir_sess)
            
            os.makedirs(dir_sess)
            
        except OSError:
            print "Warning OS error, dir_sess %s already exists"%(dir_sess)

        ##### compute envelop
        
        sess_bip_ampl_sigs,sess_timings = return_bip_ampl_sigs(dir_sess,run,bip_rc_ids,f_start,f_stop)
        
        list_bip_ampl_sigs.append(sess_bip_ampl_sigs)
        
        list_sess_timings.append(sess_timings)
        
        sess_bip_non_art,sess_keep_all_bip_non_art = return_bip_non_art(dir_sess,run,bip_rc_ids,sess_timings)
        
        #print np.sum(sess_non_art == 0,axis = 1)
        
        list_bip_non_art.append(sess_bip_non_art)
        
        list_keep_all_bip_non_art.append(sess_keep_all_bip_non_art)

    subj_bip_ampl_sigs = np.concatenate(tuple(list_bip_ampl_sigs),axis = 1)
        
    #print subj_ampl_sigs.shape
    
    #subj_sess_timings = np.concatenate(list_sess_timings,axis = 0)
    
    subj_bip_non_art = np.concatenate(tuple(list_bip_non_art),axis = 1)
    
    #print subj_non_art.shape
    
    
    
    ############# ajout keep par nom des paires d'electrodes #########################
    
    keep_pair_bip = np.zeros(shape = (len(bip_rc_ids)), dtype = bool)
    
    #print keep_pair_bip
    
    for i,bip_rc_name in enumerate(bip_rc_names):
    
        a,b = bip_rc_name.split("-")
        
        keep_pair_bip[i] = (a[0] == b[0])
        
    print np.sum(keep_pair_bip == False)
    
    
    #### filtering electrode coords and names according to number of detected artefacts
    sum_subj_bip_non_art = np.sum(subj_bip_non_art == 0,axis = 1)
    
    keep_bip_non_arts = sum_subj_bip_non_art < thr_nb_artefact
    
    print np.sum(keep_bip_non_arts == False)
    
    ########## adding to the filter rc with non registered coordinates
    
    keep_non_nans = np.sum(np.isnan(bip_rc_coords) == True,axis = 1) == 0
    
    print np.sum(keep_non_nans == False)
    
    keep_recording = np.logical_and.reduce(np.array((keep_non_nans,keep_bip_non_arts,keep_pair_bip)))
    
    print np.sum(keep_recording == False)
    
    
    ### correc_channel_names
    np_bip_rc_names = np.array(bip_rc_names,dtype = 'str')
    correc_bip_rc_names = np_bip_rc_names[keep_recording]
    print np_bip_rc_names.shape,correc_bip_rc_names.shape
    correc_bip_rc_names_file = os.path.join(dir_subj,'correct_bip_channel_names.txt')
    np.savetxt(correc_bip_rc_names_file,correc_bip_rc_names,fmt = "%s")
    
    ### correc_channel_coords
    np_bip_rc_coords = np.array(bip_rc_coords,dtype = 'float')
    print np_bip_rc_coords.shape
    correc_bip_channel_coords = np_bip_rc_coords[keep_recording,:]
    print correc_bip_channel_coords.shape
    correc_bip_rc_coords_file = os.path.join(dir_subj,'correct_bip_channel_coords.txt')
    np.savetxt(correc_bip_rc_coords_file,correc_bip_channel_coords,fmt = "%2.2f %2.2f %2.2f")
    
    
    ### logging
    with open(log_file,'w') as log:
        log.write("Removed plot for artefacts:\n")
        np.savetxt(log,np.array(zip(np_bip_rc_names[keep_bip_non_arts == False],map(str,sum_subj_bip_non_art[keep_bip_non_arts == False])),dtype = 'str'),fmt = "%s")
        
        log.write("Removed plot for void coords:\n")
        np.savetxt(log,np_bip_rc_names[keep_non_nans == False],fmt = "%s")
        
        log.write("Removed plot for non congruent electrode:\n")
        np.savetxt(log,np_bip_rc_names[keep_pair_bip == False],fmt = "%s")
        log.write("\n")
        
    #### correct_full_ampl and 
    for i,sess_index in enumerate(sess_indexes):
        
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
        correc_rc_coords_file = os.path.join(dir_sess,'correct_channel_coords.txt')
        np.savetxt(correc_rc_coords_file,correc_bip_channel_coords,fmt = "%2.2f %2.2f %2.2f")
        
        #### correct_full_ampl
        correct_bip_ampl_sigs_file = os.path.join(dir_sess,'np_correct_bip_ampl_sigs.npy')
        if not os.path.exists(correct_bip_ampl_sigs_file):
            sess_full_bip_ampl_sigs = list_bip_ampl_sigs[i][keep_recording,:][:,list_keep_all_bip_non_art[i]]
            print sess_full_bip_ampl_sigs.shape,list_bip_ampl_sigs[i].shape
            np.save(correct_bip_ampl_sigs_file,sess_full_bip_ampl_sigs)
        
        #### ampl_trigs
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
        sess_correct_trigs,sess_correct_rest,nb_removed_trigs_by_trial = return_sess_correct_trigs_and_rest(run,list_sess_timings[i],list_keep_all_bip_non_art[i],trial_types)
        with open(log_file,'a') as log:
            log.write("Session %s:\n"%sess_index)
            log.write("Number of removed artefact trigs: \n%s\n\n"%(" ".join([": ".join(trial) for trial in zip(trial_types,map(str,nb_removed_trigs_by_trial))])))
        print np.array(sess_correct_trigs).shape
        return_sess_correct_ampl_trigs(dir_sess,sess_correct_trigs,list_bip_ampl_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
        return_sess_correct_ampl_rest(dir_sess,sess_correct_rest,list_bip_ampl_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
        
def compute_subj_correct_ampl(export_dir, f_start, f_stop, name = 'SEMC',exp = 'R',thr_nb_artefact = 120):

    subject = session.query(Subject).filter_by(name = name).first()
    
    rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    rcs = [ ]
    for rcg in rcgs:
        rcs.extend(rcg.recordingchannels)
    
    
    
    dir_subj = os.path.join(export_dir,name)
    
    log_file = os.path.join(dir_subj,'log.txt')
    
    try :
        print "Creating dir_subj %s"%(dir_subj)
        
        os.makedirs(dir_subj)
        
    except OSError:
        print "Warning OS error, dir_subj %s already exists"%(dir_subj)


    #### electrode names and electrodes
    rc_names,rc_coords,rc_descriptions = return_electrode_names_and_coords(rcs,dir_subj)

    #### running
    list_ampl_sigs = []
    
    list_sess_timings = []
    
    list_keep_all_non_art = []
    
    list_non_art = []
    
    if exp == 'R':
    
        sess_indexes = ['1','2','3']
        
        trial_types = trial_R_conds
        
    elif exp == 'E':
        
        sess_indexes = ['2']
        
        trial_types = trial_E_conds
        
    
    for sess_index in sess_indexes:
        
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
            
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
        
        try :
            print "Creating dir_sess %s"%(dir_sess)
            
            os.makedirs(dir_sess)
            
        except OSError:
            print "Warning OS error, dir_sess %s already exists"%(dir_sess)

        ##### compute envelop
        
        sess_ampl_sigs,sess_timings = return_ampl_sigs(dir_sess,run,rcs,f_start,f_stop)
        
        list_ampl_sigs.append(sess_ampl_sigs)
        
        list_sess_timings.append(sess_timings)
        
        sess_non_art,sess_keep_all_non_art = return_non_art(dir_sess,run,rcs,sess_timings)
        
        #print np.sum(sess_non_art == 0,axis = 1)
        
        list_non_art.append(sess_non_art)
        
        list_keep_all_non_art.append(sess_keep_all_non_art)

    subj_ampl_sigs = np.concatenate(tuple(list_ampl_sigs),axis = 1)
        
    #print subj_ampl_sigs.shape
    
    #subj_sess_timings = np.concatenate(list_sess_timings,axis = 0)
    
    subj_non_art = np.concatenate(tuple(list_non_art),axis = 1)
    
    #print subj_non_art.shape
    
    #### filtering electrode coords and names according to number of detected artefacts
    sum_subj_non_art = np.sum(subj_non_art == 0,axis = 1)
    
    print np.sort(sum_subj_non_art)
    
    
    
    keep_non_arts = sum_subj_non_art < thr_nb_artefact
    
    print keep_non_arts
    
    ########## adding to the filter rc with non registered coordinates
    
    keep_non_nans = np.sum(np.isnan(rc_coords) == True,axis = 1) == 0
    
    print keep_non_nans
    
    keep_recording = np.logical_and(keep_non_nans,keep_non_arts)
    
    print keep_recording
    
    ### correc_channel_names
    np_rc_names = np.array(rc_names,dtype = 'str')
    
    correc_rc_names = np_rc_names[keep_recording]
    
    print np_rc_names.shape,correc_rc_names.shape
    
    
    
    
    with open(log_file,'w') as log:
        
        log.write("Removed electrodes for artefacts:\n")
        
        np.savetxt(log,np.array(zip(np_rc_names[keep_non_arts == False],map(str,sum_subj_non_art[keep_non_arts == False])),dtype = 'str'),fmt = "%s")
    
        log.write("Removed electrodes for void coords:\n")
        
        np.savetxt(log,np_rc_names[keep_non_nans == False],fmt = "%s")
    
        log.write("\n")
        
    correc_rc_names_file = os.path.join(dir_subj,'correct_channel_names.txt')
    
    np.savetxt(correc_rc_names_file,correc_rc_names,fmt = "%s")
    
    ### correc_channel_coords
    
    np_rc_coords = np.array(rc_coords,dtype = 'float')
    
    print np_rc_coords.shape
    
    correc_channel_coords = np_rc_coords[keep_recording,:]
    
    print correc_channel_coords.shape
    
    correc_rc_coords_file = os.path.join(dir_subj,'correct_channel_coords.txt')
    
    np.savetxt(correc_rc_coords_file,correc_channel_coords,fmt = "%2.2f %2.2f %2.2f")
    
    ### correc_rc_descriptions
    
    correc_rc_descriptions_file = os.path.join(dir_subj,'correct_channel_descriptions.txt')
    
    np_rc_descriptions = np.array(rc_descriptions,dtype = 'str')
    
    correc_rc_descriptions = np_rc_descriptions[keep_recording]
    
    np.savetxt(correc_rc_descriptions_file,correc_rc_descriptions,fmt = "%s")
    
    
    
    
    #### correct_full_ampl and 
    for i,sess_index in enumerate(sess_indexes):
        
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
            
        correc_rc_coords_file = os.path.join(dir_sess,'correct_channel_coords.txt')
        
        np.savetxt(correc_rc_coords_file,correc_channel_coords,fmt = "%2.2f %2.2f %2.2f")
        
        #### correct_full_ampl
     
        correct_ampl_sigs_file = os.path.join(dir_sess,'np_correct_ampl_sigs.npy')
    
        if not os.path.exists(correct_ampl_sigs_file):
        
            sess_full_ampl_sigs = list_ampl_sigs[i][keep_recording,:][:,list_keep_all_non_art[i]]
            
            print sess_full_ampl_sigs.shape,list_ampl_sigs[i].shape
            
            np.save(correct_ampl_sigs_file,sess_full_ampl_sigs)
        
            
        #### ampl_trigs
        
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
           
        sess_correct_trigs,sess_correct_rest,nb_removed_trigs_by_trial = return_sess_correct_trigs_and_rest(run,list_sess_timings[i],list_keep_all_non_art[i],trial_types)
        
        
        
        with open(log_file,'a') as log:
            
            log.write("Session %s:\n"%sess_index)
            
            log.write("Number of removed artefact trigs: \n%s\n\n"%(" ".join([": ".join(trial) for trial in zip(trial_types,map(str,nb_removed_trigs_by_trial))])))
            
        print np.array(sess_correct_trigs).shape
        
        return_sess_correct_ampl_trigs(dir_sess,sess_correct_trigs,list_ampl_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
        
        return_sess_correct_ampl_rest(dir_sess,sess_correct_rest,list_ampl_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
        
    
def compute_subj_correct_proba_density(export_dir, f_start, f_stop, name = 'SEMC',exp = 'R',thr_nb_artefact = 120):

    subject = session.query(Subject).filter_by(name = name).first()
    
    rcgs = [ rcg for rcg in subject.blocks[0].recordingchannelgroups if rcg.name.startswith('Group') ]
    rcs = [ ]
    for rcg in rcgs:
        rcs.extend(rcg.recordingchannels)
    
    
    
    dir_subj = os.path.join(export_dir,name)
    
    log_file = os.path.join(dir_subj,'log.txt')
    
    
    if exp == 'R':
    
        sess_indexes = ['1','2','3']
        
        trial_types = trial_R_conds
        
    elif exp == 'E':
        
        sess_indexes = ['2']
        
        trial_types = trial_E_conds
        
    
    
    try :
        print "Creating dir_subj %s"%(dir_subj)
        
        os.makedirs(dir_subj)
        
    except OSError:
        print "Warning OS error, dir_subj %s already exists"%(dir_subj)


    #### electrode names and electrodes
    rc_names,rc_coords,rc_descriptions = return_electrode_names_and_coords(rcs,dir_subj)

    #### running
    list_neo_sigs = []
    
    list_sess_timings = []
    
    list_keep_all_non_art = []
    
    list_non_art = []
    
    for sess_index in sess_indexes:
        
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
            
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
        
        try :
            print "Creating dir_sess %s"%(dir_sess)
            
            os.makedirs(dir_sess)
            
        except OSError:
            print "Warning OS error, dir_sess %s already exists"%(dir_sess)

        ##### compute envelop
        
        sess_neo_sigs,sess_timings = compute_neo_sigs(run,rcs)
        
        #print sess_neo_sigs
        
        list_neo_sigs.append(sess_neo_sigs)
        
        list_sess_timings.append(sess_timings)
        
        sess_non_art,sess_keep_all_non_art = return_non_art(dir_sess,run,rcs,sess_timings)
        
        #print np.sum(sess_non_art == 0,axis = 1)
        
        list_non_art.append(sess_non_art)
        
        list_keep_all_non_art.append(sess_keep_all_non_art)

    #subj_neo_sigs = np.concatenate(tuple(list_neo_sigs),axis = 1)
        
    #print subj_neo_sigs.shape
    
    
    #subj_sess_timings = np.concatenate(list_sess_timings,axis = 0)
    
    subj_non_art = np.concatenate(tuple(list_non_art),axis = 1)
    
    #print subj_non_art.shape
    
    #### filtering electrode coords and names according to number of detected artefacts
    sum_subj_non_art = np.sum(subj_non_art == 0,axis = 1)
    
    keep_non_arts = sum_subj_non_art < thr_nb_artefact
    
    print keep_non_arts
    
    ########## adding to the filter rc with non registered coordinates
    
    keep_non_nans = np.sum(np.isnan(rc_coords) == True,axis = 1) == 0
    
    print keep_non_nans
    
    keep_recording = np.logical_and(keep_non_nans,keep_non_arts)
    
    print keep_recording
    
    ### correc_channel_names
    np_rc_names = np.array(rc_names,dtype = 'str')
    
    correc_rc_names = np_rc_names[keep_recording]
    
    print np_rc_names.shape,correc_rc_names.shape
    
    with open(log_file,'w') as log:
        
        log.write("Removed electrodes for artefacts:\n")
        
        np.savetxt(log,np.array(zip(np_rc_names[keep_non_arts == False],map(str,sum_subj_non_art[keep_non_arts == False])),dtype = 'str'),fmt = "%s")
    
        log.write("Removed electrodes for void coords:\n")
        
        np.savetxt(log,np_rc_names[keep_non_nans == False],fmt = "%s")
    
        log.write("\n")
        
    correc_rc_names_file = os.path.join(dir_subj,'correct_channel_names.txt')
    
    np.savetxt(correc_rc_names_file,correc_rc_names,fmt = "%s")
    
    ### correc_channel_coords
    
    np_rc_coords = np.array(rc_coords,dtype = 'float')
    
    print np_rc_coords.shape
    
    
    correc_channel_coords = np_rc_coords[keep_recording,:]
    
    print correc_channel_coords.shape
    
    correc_rc_coords_file = os.path.join(dir_subj,'correct_channel_coords.txt')
    
    np.savetxt(correc_rc_coords_file,correc_channel_coords,fmt = "%2.2f %2.2f %2.2f")
    
    ### correc_rc_descriptions
    
    correc_rc_descriptions_file = os.path.join(dir_subj,'correct_channel_descriptions.txt')
    
    print rc_descriptions
    
    np_rc_descriptions = np.array(rc_descriptions,dtype = 'str')
    
    correc_rc_descriptions = np_rc_descriptions[keep_recording]
    
    np.savetxt(correc_rc_descriptions_file,correc_rc_descriptions,fmt = "%s")
    
    
    
    #### correct_full_ampl and 
    for i,sess_index in enumerate(sess_indexes):
        
        print sess_index
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
            
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
           
        sess_correct_odor,sess_correct_rest,nb_removed_trigs_by_trial = return_sess_correct_trigs_and_rest(run,list_sess_timings[i],list_keep_all_non_art[i],trial_types)
        
        sess_correct_neo_sigs = [neo_sig for i,neo_sig in enumerate(list_neo_sigs[i]) if keep_recording[i] == True]
        
        #print len(sees_correct_neo_sigs)
        #print sees_correct_neo_sigs[0]
        
        #### proba density
        print 'sess_proba_density_odor'
        return_sess_proba_density_odor(dir_sess,sess_correct_odor,sess_correct_neo_sigs,trial_types,f_start,f_stop)
        
        print 'sess_proba_density_rest'
        return_sess_proba_density_rest(dir_sess,sess_correct_rest,sess_correct_neo_sigs,trial_types,f_start,f_stop)
        
        
###################################################################################### testing #################################################################
def test_correct_ampl():
    
    from params import main_path,exp,envelop_method
    
    base_dir = os.path.join(main_path,"tmp-TS_" + exp + "_" + envelop_method + "_beta_all_cond_by_block_trigs")
    
    compute_subj_correct_ampl(base_dir,f_start = 15., f_stop = 35.)
    
def test_correct_proba_density():

    from params import main_path,exp,envelop_method
    
    f_start = 1.0
    f_stop = 150.0
    
    base_dir = os.path.join(main_path,"tmp-TS_" + exp + "_" + envelop_method + "_" + str(int(f_start)) + "-" + str(int(f_stop)) + "Hz_all_cond_by_block_trigs")
    
    compute_subj_correct_proba_density(base_dir,name = "VACJ",f_start = f_start, f_stop = f_stop)
    

def test_correct_bip_ampl():
    
    from params import main_path,exp,envelop_method
    
    base_dir = os.path.join(main_path,"tmp-bip_TS_" + exp + "_" + envelop_method + "_beta_all_cond_by_block_trigs")
    
    compute_subj_correct_bip_ampl(base_dir,f_start = 15., f_stop = 35.)

########################################################################## full computation ########################################
def compute_correct_ampl_all_bands():

    from params import main_path,exp,envelop_method
    
    from params import thr_nb_artefacts,freq_bands,freq_band_names

    for j,band in enumerate(freq_bands):
            
        print band
        
        base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_"+ freq_band_names[j] +"_all_cond_by_block_trigs")

        print base_dir
        
        #for i,subj in enumerate(subject_ids):
        for i,subj in enumerate(['SEMC','LEFC']):
        
            #compute_subj_correct_ampl(export_dir = base_dir,f_start = band[0], f_stop = band[1], name = subj, exp = exp,thr_nb_artefact = thr_nb_artefacts[i])
            
            compute_subj_correct_ampl(export_dir = base_dir,f_start = band[0], f_stop = band[1], name = subj, exp = exp,thr_nb_artefact = 200)
            
def compute_correct_bip_ampl_all_bands():

    from params import main_path,exp,envelop_method
    
    from params import thr_nb_artefacts,freq_bands,freq_band_names

    for j,band in enumerate(freq_bands):
            
        print band
        
        base_dir = os.path.join(main_path,"bip_TS_" + exp + "_" + envelop_method + "_"+ freq_band_names[j] +"_all_cond_by_block_trigs")

        print base_dir
        
        for i,subj in enumerate(subject_ids):
        
            compute_subj_correct_bip_ampl(export_dir = base_dir,f_start = band[0], f_stop = band[1], name = subj, exp = exp,thr_nb_artefact = thr_nb_artefacts[i])
            
def compute_correct_ampl():
    
    from params import thr_nb_artefacts,base_dir,f_start,f_stop
    
    for i,subj in enumerate(subject_ids):
        
        compute_subj_correct_ampl(export_dir = base_dir, f_start = f_start, f_stop = f_stop, name = subj,exp = 'R',thr_nb_artefact = thr_nb_artefacts[i])
        
def compute_correct_proba_density():

    from params import main_path,exp,envelop_method
    #from params import thr_nb_artefacts
    
    f_start = 1.0
    f_stop = 150.0
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + str(int(f_start)) + "-" + str(int(f_stop)) + "Hz_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    print base_dir
    
    thr_nb_artefacts = [150]
    #thr_nb_artefacts = [150,200]
    
    for i,subj in enumerate(['SEMC']):
    #for i,subj in enumerate(['SEMC','LEFC']):
    
        compute_subj_correct_proba_density(export_dir = base_dir,f_start = f_start, f_stop = f_stop, name = subj, exp = exp,thr_nb_artefact = thr_nb_artefacts[i])
        

if __name__ == '__main__':
    
    ### test tout d'un coup
    #test_correct_ampl()
    #test_correct_bip_ampl()
    
    #test_correct_proba_density()
    
    ### mise en place pipeline:
    #compute_correct_ampl()
    #compute_correct_ampl_all_bands()
    #compute_correct_bip_ampl_all_bands()
    
    compute_correct_proba_density()
    
    
    
    