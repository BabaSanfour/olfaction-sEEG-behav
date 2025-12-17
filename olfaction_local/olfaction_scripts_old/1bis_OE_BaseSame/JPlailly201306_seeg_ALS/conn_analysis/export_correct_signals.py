# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 
sys.path.append('../seeg_analysis') 

from connection import *
from params import *
from OpenElectrophy.timefrequency import TimeFreq

from artifact_detection import is_artifact_in_window,detect_artifact_on_one_sig,return_artifact_timings,clean_trig

timefreqtools.set_session(session, dbinfo)
  
import itertools
import pandas as pd

from matplotlib import pyplot
from export_inband_amplitude import return_electrode_names_and_coords,return_non_art,return_sess_correct_trigs_and_rest


def compute_sigs(run,rcs):
    
    sigs = []
    
    if len(rcs) != 0:
        
        for i in range(len(rcs)):
            
            print i
            
            anasig_i = run.segments[0].analogsignals.filter_by(recordingchannel_id = rcs[i].id)[0]
            
            anasig = session.query(AnalogSignal).get(anasig_i.id)
            
            neosig = anasig.to_neo()
            
            times = neosig.times.rescale('s').magnitude
            
            sig = neosig.magnitude
            
            sigs.append(sig)
            
    else:
        
        times = np.array([])
    #print sigs
    
    return np.array(sigs),times
    
    
    
#####################" save full unsegmented signals (npy)
def return_sigs(cur_dir,run,rcs):
    
    np_sigs_file = os.path.join(cur_dir,'np_sigs.npy')
    
    timings_file = os.path.join(cur_dir,'timings.npy')
        
    if not os.path.exists(np_sigs_file) or not os.path.exists(timings_file):
        
        np_sigs,timings = compute_sigs(run,rcs)
        
        np.save(np_sigs_file,np_sigs)
        
        np.save(timings_file,timings)
    
    else:
        
        np_sigs = np.load(np_sigs_file)
        
        timings = np.load(timings_file)
        
    return np_sigs,timings

#~ def return_mat_sigs(cur_dir,run,rcs):
    
    #~ mat_sigs_file = os.path.join(cur_dir,'np_sigs.mat')
    
    #~ if not os.path.exists(mat_sigs_file):
        
        #~ np_sigs,timings = compute_sigs(run,rcs)
        
        #~ savemat(mat_sigs_file,{'sigs': np_sigs,'timings': timings})
        
    #~ else:
        
        #~ mat_sigs = loadmat(mat_sigs_file)
        
        #~ np_sigs = mat_sigs['sigs']
        
        #~ timings = mat_sigs['timings']        
        
        #~ timings = timings.reshape(-1)
        
    #~ print np_sigs
    
    #~ return np_sigs,timings


############################ timings ########################################

############ compute timings (from times)
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
        
################## save timings (npy)
def return_timings(timings_file,run,rcs):

    if not os.path.exists(timings_file):
        
        timings = compute_timings(run,rcs)
        
        np.save(timings_file,timings)
                        
        
    else:
        
        timings = np.load(timings_file)
        
    return timings

    
############################# calcule la segmentation (en fonction des bons trigs) des signaux corrects apres rejection artefact 
def compute_trial_sigs(correct_sigs,correct_trigs,sess_timings):

    trial_sigs = []
    
    print correct_trigs
    
    old_len_sel = -1
    
    for trig in correct_trigs:
        
        #print trig
        
        sel = (trig[0] < sess_timings) & (sess_timings <= trig[1])
        
        if old_len_sel == -1:
            
            old_len_sel = np.sum(sel == True)
            
        print np.sum(sel == True) 
        
        if np.sum(sel == True) != old_len_sel:
            
            print trig
            
            print trig[1],round(trig[1],3)
            
            sel = (trig[0] < sess_timings) & (sess_timings <= round(trig[1],3))
        
        if np.sum(sel == True) != 0:

            print correct_sigs[:,sel].shape
            
            trial_sigs.append(correct_sigs[:,sel])
            
            old_len_sel = np.sum(sel == True)
        
        else:
            
            print "Warning, trig %f %f is out of timings, skipping"%(trig[0], trig[1]) 
    return np.array(trial_sigs)

    
    
############################# sauve les signaux (trigs) correct triggés apres rejection artefact au format npy (peut etre fait au format .mat) 
def return_sess_correct_trigs(write_path,sess_correct_trigs,sess_correct_sigs,sess_timings,trial_types):

    print len(sess_correct_sigs)
    
    sess_correct_sig_trigs = []
    
    for i,trial_type in enumerate(trial_types):
        
        #### npy
        np_trigs_file = os.path.join(write_path,'correct_ts_' + str(trial_type) + '_by_odor_trigs.npy')
        
        if not os.path.exists(np_trigs_file):
            
            trial_sigs = compute_trial_sigs(sess_correct_sigs,sess_correct_trigs[i],sess_timings)
            
            print trial_sigs.shape
            
            np.save(np_trigs_file,trial_sigs)
        else:
            
            trial_sigs = np.load(np_trigs_file)
            
        ### mat
        #~ mat_trigs_file = os.path.join(write_path,'correct_ts_' + trial_type + '_by_odor_trigs.mat')
        
        #~ if not os.path.exists(mat_trigs_file):
            
            #~ trial_sigs = compute_trial_sigs(sess_correct_sigs,sess_correct_trigs[i],sess_timings)
            
            #~ print trial_sigs.shape
            
            #~ savemat(mat_trigs_file,{'sigs':trial_sigs})
        #~ else:
            
            #~ trial_sigs = loadmat(mat_trigs_file)['sigs']
           
           
        sess_correct_sig_trigs.append(trial_sigs)
        
    return sess_correct_sig_trigs
            
############################# sauve les signaux (rest) correct triggés apres rejection artefact au format npy (peut etre fait au format .mat) 
def return_sess_correct_rest(write_path,sess_correct_trigs,sess_correct_sigs,sess_timings,trial_types):

    print len(sess_correct_sigs)
    
    sess_correct_sig_rest = []
    
    for i,trial_type in enumerate(trial_types):
        
        ####npy
        rest_file = os.path.join(write_path,'correct_ts_' + str(trial_type) + '_by_rest_trigs.npy')
        
        if not os.path.exists(rest_file):
            
            trial_sigs = compute_trial_sigs(sess_correct_sigs,sess_correct_trigs[i],sess_timings)
            
            print trial_sigs.shape
            
            np.save(rest_file,trial_sigs)
        else:
            
            trial_sigs = np.load(rest_file)
        
        ### mat
        #~ rest_file = os.path.join(write_path,'correct_ts_' + trial_type + '_by_rest_trigs.mat')
        
        #~ if not os.path.exists(rest_file):
            
            #~ trial_sigs = compute_trial_sigs(sess_correct_sigs,sess_correct_trigs[i],sess_timings)
            
            #~ print trial_sigs.shape
            
            #~ savemat(rest_file,{'sigs':trial_sigs})
        #~ else:
            
            #~ trial_sigs = loadmat(rest_file)['sigs']
            
        sess_correct_sig_rest.append(trial_sigs)
        
    return sess_correct_sig_rest
    
############################# selection des signaux pour toute la session
def return_sess_correct_win(dir_sess,run,sess_timings,sess_keep_all_non_art,trial_types):


    #### compute all correct trigs and rest periods 
    sess_correct_win = []
    
    sess_correct_trial_ids = []
    
    nb_removed_trigs_by_trial = []
    
    for trial_type in trial_types:
    #for trial_type in ['hit']:
        
        print trial_type
        
        ## odor + rest
        correct_win,correct_trial_ids,nb_removed_trigs = oect_win(run,sess_timings,sess_keep_all_non_art,trial_type = trial_type)
            
        print correct_trial_ids

        
        correct_trial_ids_file = os.path.join(dir_sess, "correct_trial_ids_win_"+ trial_type + ".txt")

        with open(correct_trial_ids_file,'w') as f:
            json.dump(correct_trial_ids,f)

        sess_correct_win.append(correct_win)
        
        sess_correct_trial_ids.append(correct_trial_ids)
        
        nb_removed_trigs_by_trial.append(nb_removed_trigs)
        
    #print sess_correct_trigs
    
    return sess_correct_win,sess_correct_trial_ids,nb_removed_trigs_by_trial
    
    
######################################################### main analyses ##################################################################################################################################


def compute_subj_correct_sigs(export_dir, name = 'SEMC',exp = 'E',thr_nb_artefact = 120):

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


    if exp == 'R':
    
        sess_indexes = ['1','2','3']
        
        trial_types = trial_R_conds
        
    elif exp == 'E':
        
        sess_indexes = ['1','2']
        
        trial_types = trial_E_conds
        
        
    list_sigs = []
    
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
        
        ### npy
        sess_sigs,sess_timings = return_sigs(dir_sess,run,rcs)
        
        ### mat
        #sess_sigs,sess_timings = return_mat_sigs(dir_sess,run,rcs)
        
        
        #sess_ampl_sigs,sess_timings = return_ampl_sigs(dir_sess,run,rcs,f_start,f_stop)
        
        print sess_sigs.shape
        
        print sess_timings.shape
        
        list_sigs.append(sess_sigs)
        
        list_sess_timings.append(sess_timings)
 
        sess_non_art,sess_keep_all_non_art = return_non_art(dir_sess,run,rcs,sess_timings)
        
        #print np.sum(sess_non_art == 0,axis = 1)
        
        list_non_art.append(sess_non_art)
        
        list_keep_all_non_art.append(sess_keep_all_non_art)


    subj_sigs = np.concatenate(tuple(list_sigs),axis = 1)
        
    subj_non_art = np.concatenate(tuple(list_non_art),axis = 1)
    
    #print subj_non_art.shape
    
    #### filtering electrode coords and names according to number of detected artefacts
    sum_subj_non_art = np.sum(subj_non_art == 0,axis = 1)
    print sum_subj_non_art
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
    
    ### correc_all_info
    #path = r'C:/Users/Anne-Lise/Dropbox/Intra_EM/1bis_OE_BaseSam/TS_R_all_cond_by_block_trigs/CHAF/'
    path = os.path.join(dir_subj)
    x = np.loadtxt(path + '/' + 'correct_channel_names.txt', dtype=str)
    names_frame = pd.DataFrame(x)
    y = np.loadtxt(path + '/' +'correct_channel_coords.txt', dtype=str)
    coord_frame = pd.DataFrame(y)
    z = np.loadtxt(path + '/' +'correct_channel_descriptions.txt', dtype=str)
    descr_frame = pd.DataFrame(z)
    all_df = [names_frame,coord_frame,descr_frame ]
    all_info = pd.concat(all_df, axis=1)
    print all_info, all_info.shape
    #correc_all_info_file = os.path.join(dir_subj,'correct_all_info.txt')
    all_info.to_csv(path+"correct_all_info.txt", header=False, index=False)
         
    
    #### correct_full_ampl and 
    for i,sess_index in enumerate(sess_indexes):
        
        dir_sess = os.path.join(dir_subj , exp + str(sess_index))
            
        correc_rc_coords_file = os.path.join(dir_sess,'correct_channel_coords.txt')
        
        np.savetxt(correc_rc_coords_file,correc_channel_coords,fmt = "%2.2f %2.2f %2.2f")
        
        #### correct_full_ampl
        ### npy
        
        correct_sigs_file = os.path.join(dir_sess,'np_correct_sigs.npy')
    
        if not os.path.exists(correct_sigs_file):
        
            #sess_full_sigs = list_sigs[i][keep_recording,:][:,list_keep_all_non_art[i]]
            sess_full_sigs = list_sigs[i]
            print sess_full_sigs.shape,list_sigs[i].shape
            
            np.save(correct_sigs_file,sess_full_sigs)
        
        ### mat
        #~ mat_correct_sigs_file = os.path.join(dir_sess,'correct_sigs.mat')
    
        #~ if not os.path.exists(mat_correct_sigs_file):
        
            #~ sess_full_sigs = list_sigs[i][keep_recording,:][:,list_keep_all_non_art[i]]
            
            #~ print sess_full_sigs.shape,list_sigs[i].shape
            
            #~ savemat(mat_correct_sigs_file,{'sigs':sess_full_sigs})
        
            
        #### ampl_trigs
        
        run = subject.runs.filter_by(exp= exp).filter_by(index=sess_index).first()
           
        sess_correct_trigs,sess_correct_rest,nb_removed_trigs_by_trial = return_sess_correct_trigs_and_rest(run,list_sess_timings[i],list_keep_all_non_art[i],trial_types)
        
        
        
        with open(log_file,'a') as log:
            
            log.write("Session %s:\n"%sess_index)
            
            log.write("Number of removed artefact trigs: \n%s\n\n"%(" ".join([": ".join(trial) for trial in zip(str(trial_types),map(str,nb_removed_trigs_by_trial))])))
            
        print np.array(sess_correct_trigs).shape
        
        return_sess_correct_trigs(dir_sess,sess_correct_trigs,list_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
        
        return_sess_correct_rest(dir_sess,sess_correct_rest,list_sigs[i][keep_recording,:],list_sess_timings[i],trial_types)
                
##################################################### testing #################################################################

def test_correct_sigs():
    
    from params import main_path,exp
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_all_cond_by_block_trigs")
    
    
    compute_subj_correct_sigs(base_dir, name = 'SEMC', exp = exp)
    

################################################## full computation #######################################
def compute_all_correct_sigs():

    from params import main_path,exp
    
    from params import thr_nb_artefacts

    for exp in ['E']:
            
        base_dir = os.path.join(main_path,"TS_"+exp+"_all_cond_by_block_trigs")

        print base_dir
        
        for i,subj in enumerate(subject_ids):
        #for i,subj in enumerate(['PIRJ']):
        
            compute_subj_correct_sigs(base_dir,name = subj, exp = exp,thr_nb_artefact = thr_nb_artefacts[i])
            
if __name__ == '__main__':
    
    ### test
    #test_correct_sigs()
    
    ### mise en place pipeline:
    compute_all_correct_sigs()
    
    
    
