# -*- coding: utf-8 -*-
import sys, os
sys.path.append('../behavior') 
from connection import *
import getpass

database_name = url.split('/')[-1]
if sys.platform.startswith('win')and getpass.getuser() == 'alsaive':
    #~ desktop = 'C:/Users/{}/Desktop'.format(getpass.getuser())
    main_path = "D:/Episodic_Epilepsy/Data-seeg"
    
elif sys.platform.startswith('linux') and getpass.getuser() in [ 'sgarcia', 'samuel' ]:
    desktop = '/home/{}/Bureau'.format(getpass.getuser())
    main_path = desktop

elif sys.platform.startswith('linux') and getpass.getuser() == 'david':
    desktop = '/home/{}/Bureau'.format(getpass.getuser())
    main_path = "/mnt/Results/AnneLise-Intra-episodic"
    data_path = "/mnt/Results/AnneLise-Intra-episodic"




#if not os.path.exists(joblib_cachedir):
    #os.mkdir(joblib_cachedir)  
    
sfreq = 512

    
baseline_mode = 'no-normalisation'
#~ baseline_mode = 'ratio'
#~ baseline_mode = 'difference'
#~ baseline_mode = 'z-score'


#envelop_method = 'hilbert'
envelop_method = 'timefreq'

#t_win_start_odor = -0.5
#t_win_start_odor = 0.
t_win_start_odor = 0.5
#t_win_start_odor = 1.0

#t_win_stop_odor = 0.0
#t_win_stop_odor = 0.5
t_win_stop_odor = 1.0
#t_win_stop_odor = 1.5

t_win_start_rest = 0.5
t_win_stop_rest = 1.0
#t_win_stop_rest = 1.5

subject_ids = ['FERJ','SEMC','CHAF','LEFC','VACJ']
thr_nb_artefacts = [200,150,200,200,200,]

#subject_ids = ['VACJ']
#thr_nb_artefacts = [200]

#subject_ids = ['SEMC']
#thr_nb_artefacts = [150]

exp = 'E'
#exp = 'R'



#trial_R_conds = ['hit','cr','fa','miss','all']
#trial_R_conds = ['hit','cr','fa','miss']
trial_R_conds = ['all']
trial_E_conds = ['all']



if exp == 'R':

    sess_indexes = ['1','2','3']
    
    trial_types = trial_R_conds
    
elif exp == 'E':
    
    sess_indexes = ['1','2']
    
    trial_types = trial_E_conds
    


### "gamma" band
#f_start = 60.
#f_stop = 90.

#### "beta" band
#f_start = 15.
#f_stop = 35.

### "delta"
#f_start = 4.
#f_stop = 8.

#base_dir = os.path.join(main_path,"TS_R_all_cond")

#base_dir = os.path.join(main_path,"TS_R_" + envelop_method +"_all_cond_by_trial")

#base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_"+str(int(f_start)) + "-"+ str(int(f_stop)) +"_all_cond_by_block_trigs")

#correl_rest_analysis_name = "correl_inband_ampl_by_rest_trial_" + envelop_method
#correl_analysis_name = "correl_inband_ampl_by_odor_trial_" + envelop_method

#correl_analysis_name = "correl_inband_ampl_by_odor_trigs_"+str(int(t_win_start)) + "-" +str(int(t_win_stop)) +"s_" + envelop_method
#correl_rest_analysis_name = "correl_inband_ampl_by_rest_trigs_"+str(int(t_win_start)) + "-" +str(int(t_win_stop)) +"s_" + envelop_method

#correl_analysis_name = "correl_inband_ampl_by_block_odor_trigs_"+str(int(t_win_start)) + "-" +str(int(t_win_stop)) +"s_" + envelop_method + "_"+str(int(f_start)) + "-"+ str(int(f_stop)) 
#correl_rest_analysis_name = "correl_inband_ampl_by_rest_trigs_"+str(int(t_win_start)) + "-" +str(int(t_win_stop)) +"s_" + envelop_method + "_"+str(int(f_start)) + "-"+ str(int(f_stop)) 

### by_band

freq_bands = [[4.,8.],[15.,35.],[60.,90.]]

freq_band_names = ["theta","beta","gamma"]




con_method = "coh"

#column_label = "choosen_port_name"

column_label = "sscore_WhatWhere"


correl_analysis_name = con_method + "_sigs_" + exp + "_sep_by_trigs_by_band"



############## proba_density

##KS_alpha = 1.36 ### for p = 0.05
#KS_alpha = 1.95 ### for p = 0.001

############## correlation amplitude

### Z_thr for correlations
##Z_thr = 1.96 ### p = 0.05
#Z_thr = 2.58 ### p = 0.01
##Z_thr = 2.8 ### p = 0.005

#conf_interval_prob = 0.05

##### Radatools directory
#radatools_path ='/home/david/Packages/radatools/radatools-3.2-linux32/'
#radatools_optim = "WS trfr 100"

#### full sig
##correl_analysis_name = "correl_inband_full_ampl_by_band"

#### concat sig for trigs or odors
#correl_analysis_name = "correl_ampl_by_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")
##correl_analysis_name = "correl_ampl_by_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")

#### separate sig for trigs or odors
##correl_analysis_name = "correl_ampl_by_sep_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")
##correl_analysis_name = "correl_ampl_by_sep_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")
    
#coclass_analysis_name = "coclass_rada_by_sep_trigs_by_band"

#split_coclass_analysis_name = coclass_analysis_name.split('_')
