# -*- coding: utf-8 -*-
import sys, os
sys.path.append('../behavior')
from connection import *
import getpass

database_name = url.split('/')[-1]
if sys.platform.startswith('win')and getpass.getuser()== 'alsaive':
    main_path = "D:/Episodic_Epilepsy/Data-seeg"
    joblib_cachedir = os.path.join(main_path,'cache_seeg')

if sys.platform.startswith('linux')and getpass.getuser()== 'karim':
    main_path = "/media/karim/Datas4To"
    joblib_cachedir = os.path.join(main_path,'cache_seeg')

elif sys.platform.startswith('win')and getpass.getuser()== 'Anne-Lise':
    main_path = "C:/Users/Anne-Lise/Dropbox/Intra_EM/1bis_OE_BaseSam"
    joblib_cachedir = os.path.join(main_path,'cache_seeg')

elif sys.platform.startswith('linux') and getpass.getuser() in [ 'sgarcia', 'samuel' ]:
    desktop = '/home/{}/Bureau'.format(getpass.getuser())
    main_path = desktop
    joblib_cachedir = os.path.join(main_path, 'cache_calcul_seeg_'+database_name)

elif sys.platform.startswith('linux') and getpass.getuser() == 'david':
    desktop = '/home/{}/Bureau'.format(getpass.getuser())
    main_path = "/mnt/Data/episodic_ieeg"
    joblib_cachedir = os.path.join(main_path, 'cache_calcul_seeg_'+database_name)


#if not os.path.exists(joblib_cachedir):
    #os.mkdir(joblib_cachedir)


baseline_mode = 'no-normalisation'
#~ baseline_mode = 'ratio'
#~ baseline_mode = 'difference'
#~ baseline_mode = 'z-score'


#envelop_method = 'hilbert'
envelop_method = 'timefreq'

#t_win_start_odor = -0.5
#t_win_start_odor = 0.
t_win_start_odor = -5#-3#-1.5
#t_win_start_odor = 1.0

#t_win_stop_odor = 0.0
#t_win_stop_odor = 0.5
t_win_stop_odor = 3.#4#2
#t_win_stop_odor = 1.5

t_win_start_rest = 0. #-1.7
t_win_stop_rest = 1. #+0.7
#t_win_start_rest = -2.
#t_win_stop_rest = 0

subject_ids = ['PIRJ','SEMC','VACJ', 'LEFC','FERJ','CHAF']
#subject_ids = ['PIRJ',]
#subject_ids = ['PIRJ', 'VACJ']
#~ subject_ids = ['VACJ', 'PIRJ', 'LEFC']
thr_nb_artefacts = [200,200,200,200,200,200]
#thr_nb_artefacts = [400,400,400,400,400,400,400,400]
#thr_nb_artefacts = [300,300,300,300,300,300]

#subject_ids = ['SEMC']
#thr_nb_artefacts = [150]

exp = 'R'
#exp = 'R'

#trial_R_conds = ['hit','cr','fa','miss','all']
#trial_R_conds = ['hit','cr','fa','miss']
trial_R_conds = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']#['all']
#trial_E_conds = ['no_odor']
#trial_R_conds = ['no_odor']
#trial_E_conds = ['all']
trial_E_conds = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']

#trial_E_conds = range(1,21)
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

############# proba_density

#KS_alpha = 1.36 ### for p = 0.05
KS_alpha = 1.95 ### for p = 0.001

############# correlation amplitude

## Z_thr for correlations
#Z_thr = 1.96 ### p = 0.05
Z_thr = 2.58 ### p = 0.01
#Z_thr = 2.8 ### p = 0.005

conf_interval_prob = 0.05

#### Radatools directory
radatools_path ='/home/david/Packages/radatools/radatools-3.2-linux32/'
radatools_optim = "WS trfr 100"

### full sig
#correl_analysis_name = "correl_inband_full_ampl_by_band"

### concat sig for trigs or odors
correl_analysis_name = "correl_ampl_by_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")
#correl_analysis_name = "correl_ampl_by_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")

### separate sig for trigs or odors
#correl_analysis_name = "correl_ampl_by_sep_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")
#correl_analysis_name = "correl_ampl_by_sep_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_")

coclass_analysis_name = "coclass_rada_by_sep_trigs_by_band"

split_coclass_analysis_name = coclass_analysis_name.split('_')
