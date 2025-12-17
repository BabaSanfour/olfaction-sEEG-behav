# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 

from connection import *

from params import main_path,joblib_cachedir
from params import subject_ids,exp,trial_R_conds,trial_E_conds

from params import conf_interval_prob,Z_thr
from params import radatools_optim,radatools_path

from params import freq_band_names

from params import correl_analysis_name
    
from params import exp,envelop_method

from params import t_win_start_odor,t_win_stop_odor,t_win_start_rest,t_win_stop_rest
#baseline_t_win_start,
#baseline_t_win_stop,f_start,f_stop,envelop_method,baseline_mode

#### basic imports
import sys,io,os,fnmatch,shutil

import matplotlib
matplotlib.use('PS')

#### nibabel import
import nibabel as nib

##### nipype import
#from nipype import config
#config.enable_debug_mode()
import nipype
print nipype.__version__

import nipype.interfaces.io as nio

from nipype.interfaces.utility import IdentityInterface,Function
import nipype.pipeline.engine as pe

from dmgraphanalysis_nodes.nodes.correl_mat import ComputeConfCorMat,ConcatTS,SeparateTS
from dmgraphanalysis_nodes.nodes.modularity import ComputeNetList,PrepRada,CommRada,PlotIGraphModules

 
from dmgraphanalysis_nodes.utils import get_first


def gather_proba_density(subj_id = 'SEMC', freq = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat
    
    from plot_electrodes import hex_to_rgb
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + freq + "_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    labels_file = os.path.join(base_dir,subj_id,"correct_channel_names.txt")
    
    labels = [line.strip() for line in open(labels_file)]
    
    print labels
    
    ####### orig sig
    
    proba_density_odor_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_odor_trigs.npy")
    
    proba_density_odor = np.load(proba_density_odor_file)
    
    print proba_density_odor
    
    proba_density_rest_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_rest_trigs.npy")
    
    proba_density_rest = np.load(proba_density_rest_file)
    
    print proba_density_rest.shape
    
    list_diff_cumsum_proba = []
        
    for i,label in enumerate(labels):
        
        print i,label
        
        cumsum_proba_odor = np.cumsum(proba_density_odor[i,:,:],axis = 1) 
        
        cumsum_proba_rest = np.cumsum(proba_density_rest[i,:,:], axis = 1)
        
        cumsum_proba_odor_file = os.path.join(base_dir,subj_id,sess_index, "plot_cumsum_proba_odor_" + label + ".eps")
        plot_ranged_cormat(cumsum_proba_odor_file,cumsum_proba_odor,fix_full_range = [0.0,1.0])
        
        cumsum_proba_rest_file = os.path.join(base_dir,subj_id,sess_index, "plot_cumsum_proba_rest_" + label + ".eps")
        plot_ranged_cormat(cumsum_proba_rest_file,cumsum_proba_rest,fix_full_range = [0.0,1.0])
        
        
        ### diff
        diff_proba_density = cumsum_proba_odor - cumsum_proba_rest
        
        print diff_proba_density.shape
        
        
        diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_diff_proba_density_" + label + ".eps")
        plot_ranged_cormat(diff_proba_density_file,diff_proba_density,fix_full_range = [-0.5,0.5])
        
        list_diff_cumsum_proba.append(diff_proba_density)
        
    all_diff_cumsum_proba = np.array(list_diff_cumsum_proba)
    
    print all_diff_cumsum_proba.shape
    
    mean_all_diff_cumsum_proba = np.mean(all_diff_cumsum_proba, axis = 0)

    diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_diff_proba_density_mean_all_plots.eps")
    plot_ranged_cormat(diff_proba_density_file,mean_all_diff_cumsum_proba,fix_full_range = [-0.5,0.5])
    

def gather_proba_density_only_GM(subj_id = 'SEMC', freq = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat
    
    from plot_electrodes import hex_to_rgb
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + freq + "_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    labels_file = os.path.join(base_dir,subj_id,"correct_channel_names.txt")
    
    labels = [line.strip() for line in open(labels_file)]
    
    print labels
    
    descriptions_file = os.path.join(base_dir,subj_id,"correct_channel_descriptions.txt")
    
    descriptions = [line.strip() for line in open(descriptions_file)]
    
    print descriptions
    
    ####### orig sig
    
    proba_density_odor_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_odor_trigs.npy")
    
    proba_density_odor = np.load(proba_density_odor_file)
    
    print proba_density_odor
    
    proba_density_rest_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_rest_trigs.npy")
    
    proba_density_rest = np.load(proba_density_rest_file)
    
    print proba_density_rest.shape
    
    list_diff_cumsum_proba = []
        
    for i,label in enumerate(labels):
        
        #print descriptions[i]
        
        if descriptions[i] in ['WM','LCR','skull']:
            
            if descriptions[i] in ['skull']:
                print descriptions[i]
            
            continue
        
        
        #print i,label
        
        list_diff_cumsum_proba.append(np.cumsum(proba_density_odor[i,:,:],axis = 1) - np.cumsum(proba_density_rest[i,:,:], axis = 1))
        
    all_diff_cumsum_proba = np.array(list_diff_cumsum_proba)
    
    print all_diff_cumsum_proba.shape
    
    mean_all_diff_cumsum_proba = np.mean(all_diff_cumsum_proba, axis = 0)

    diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_diff_proba_density_mean_all_GM_plots.eps")
    plot_ranged_cormat(diff_proba_density_file,mean_all_diff_cumsum_proba,fix_full_range = [-0.5,0.5])
    
def gather_proba_density_by_description(subj_id = 'SEMC', freq = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat
    
    from plot_electrodes import hex_to_rgb
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + freq + "_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    labels_file = os.path.join(base_dir,subj_id,"correct_channel_names.txt")
    
    labels = [line.strip() for line in open(labels_file)]
    
    print labels
    
    descriptions_file = os.path.join(base_dir,subj_id,"correct_channel_descriptions.txt")
    
    descriptions = [line.strip() for line in open(descriptions_file)]
    
    print descriptions
    
    print np.unique(descriptions)
     
    ####### orig sig
    
    proba_density_odor_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_odor_trigs.npy")
    
    proba_density_odor = np.load(proba_density_odor_file)
    
    print proba_density_odor
    
    proba_density_rest_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_rest_trigs.npy")
    
    proba_density_rest = np.load(proba_density_rest_file)
    
    print proba_density_rest.shape
    
    
    log_file = os.path.join(base_dir,subj_id,sess_index,'nb_plots_by_desc.txt')
    
    nb_plots_by_desc = []
    
    for desc in np.unique(descriptions):
        
        print desc 
        
        def check_description_name(desc):
            
            split_desc = desc.split('-')
            
            print split_desc
            
            for desc2 in split_desc:
                
                if desc2 in ['WM','LCR','skull']:
                    
                    return False
            return True
            
        
        if check_description_name(desc) == False:
            
            continue
        
                
        
        list_diff_cumsum_proba = []
            
            
        pos_des, = np.where(np.array(descriptions,dtype = 'str') == desc)
        
        print pos_des
        
        nb_plots_by_desc.append((desc,str(len(pos_des))))
        
        
        
        all_diff_cumsum_proba = np.array([np.cumsum(proba_density_odor[i,:,:],axis = 1) - np.cumsum(proba_density_rest[i,:,:], axis = 1) for i in pos_des])
        
        print all_diff_cumsum_proba.shape
        
        mean_all_diff_cumsum_proba = np.mean(all_diff_cumsum_proba, axis = 0)

        diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_diff_proba_density_mean_all_" + desc + "_plots.eps")
        plot_ranged_cormat(diff_proba_density_file,mean_all_diff_cumsum_proba,fix_full_range = [-0.5,0.5])
        
    np_nb_plots_by_desc = np.array(nb_plots_by_desc, dtype = 'str')
    
    print np_nb_plots_by_desc
    
    np.savetxt(log_file,np_nb_plots_by_desc, fmt = "%s")
    
def gather_proba_density_by_chunks(subj_id = 'SEMC', freq = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from params import KS_alpha
    
    from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat
    
    from plot_electrodes import hex_to_rgb
    
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + freq + "_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    labels_file = os.path.join(base_dir,subj_id,"correct_channel_names.txt")
    
    labels = [line.strip() for line in open(labels_file)]
    
    print labels
    
    ####### orig sig
    
    proba_density_odor_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_chunks_" + trial_cond + "_by_odor_trigs.npy")
    
    proba_density_odor = np.load(proba_density_odor_file)
    
    print proba_density_odor.shape
    
    proba_density_rest_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_chunks_" + trial_cond + "_by_rest_trigs.npy")
    
    proba_density_rest = np.load(proba_density_rest_file)
    
    print proba_density_rest.shape
    
    for i,label in enumerate(labels):
        
        print i,label
        
        count_signif_diff_odor_rest = np.zeros(shape = proba_density_odor.shape[2:],dtype = 'float') 
        
        print count_signif_diff_odor_rest.shape
        
        
        ### vrai KS
        #count_signif_diff_odor_rest = np.zeros(shape = proba_density_odor.shape[2],dtype = 'int64') 
        
        #print count_signif_diff_odor_rest.shape
        
        
        
        
        list_diff_cumsum_proba = []
        
        for j in range(proba_density_odor.shape[1]):
            
            
            cumsum_proba_odor = np.cumsum(proba_density_odor[i,j,:,:],axis = 1) 
            
            cumsum_proba_rest = np.cumsum(proba_density_rest[i,j,:,:], axis = 1)
            
            
            ### diff
            diff_proba_density = cumsum_proba_odor - cumsum_proba_rest
            
            list_diff_cumsum_proba.append(diff_proba_density)
            
            print diff_proba_density.shape
            
            #### D_val (vrai KS)
            #max_diff_proba_density_by_freq = np.amax(diff_proba_density, axis = 1)
            
            #print max_diff_proba_density_by_freq.shape
            
            #KS_tresh = KS_alpha * np.sqrt(2.0/diff_proba_density.shape[1])
            
            #print KS_tresh
            
            #print KS_tresh < np.abs(max_diff_proba_density_by_freq)
            
            ### D_val
            KS_tresh = KS_alpha * np.sqrt(2.0/diff_proba_density.shape[1])
            
            print KS_tresh
            
            print KS_tresh < np.abs(diff_proba_density)
            
            count_signif_diff_odor_rest[KS_tresh < np.abs(diff_proba_density)] = count_signif_diff_odor_rest[KS_tresh < np.abs(diff_proba_density)] + 1.0 
            
            
        mean_diff_proba_density = np.mean(np.array(list_diff_cumsum_proba),axis = 0)
        
        diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_mean_chunks_diff_proba_density_" + label + ".eps")
        plot_ranged_cormat(diff_proba_density_file,diff_proba_density,fix_full_range = [-0.5,0.5])
        
        count_signif_diff_odor_rest = count_signif_diff_odor_rest/float(proba_density_odor.shape[1])
        
        print count_signif_diff_odor_rest
        
        diff_proba_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_count_signif_diff_proba_density_" + label + ".eps")
        plot_ranged_cormat(diff_proba_density_file,count_signif_diff_odor_rest,fix_full_range = [0,1.0])
        
def sum_modular_proba_density(subj_id = 'SEMC', freq = 'beta', freq_modules = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from dmgraphanalysis_nodes.utils_net import read_lol_file
    from dmgraphanalysis_nodes.utils_plot import plot_cormat
    from dmgraphanalysis_nodes.utils_igraph import igraph_colors
    
    from plot_electrodes import hex_to_rgb
    
    
    path = os.path.join(main_path, correl_analysis_name, "_freq_band_" + freq_modules + "_sess_index_" + sess_index + "_subject_id_" + subj_id + "_trial_cond_" + trial_cond)
    
    lol_file = os.path.join(path,"community_rada","Z_List.lol")
    
    community_vect = read_lol_file(lol_file)
    
    print community_vect
    
    colors_by_mod = [igraph_colors[i] for i in np.unique(community_vect)]
    
    print colors_by_mod
    
    colors_plots_by_mod = np.array(colors_by_mod,dtype = 'str')[community_vect]
    
    print colors_plots_by_mod
        
    base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_" + freq + "_all_cond_by_block_trigs_odor" + str(t_win_start_odor).replace('.','_') + "-" + str(t_win_stop_odor).replace('.','_') + "s_rest" + str(t_win_start_rest).replace('.','_') + "-" + str(t_win_stop_rest).replace('.','_') + "s")
    
    labels_file = os.path.join(base_dir,subj_id,"correct_channel_names.txt")
    
    labels = [line.strip() for line in open(labels_file)]
    
    print labels
    
    ####### orig sig
    
    proba_density_odor_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_odor_trigs.npy")
    
    proba_density_odor = np.load(proba_density_odor_file)
    
    print proba_density_odor
    
    proba_density_rest_file = os.path.join(base_dir,subj_id,sess_index,"proba_density_" + trial_cond + "_by_rest_trigs.npy")
    
    proba_density_rest = np.load(proba_density_rest_file)
    
    print proba_density_rest.shape
    
    
    list_diff_cumsum_proba = []
        
    for i,label in enumerate(labels):
        
        print i,label
        
        cumsum_proba_odor = np.cumsum(proba_density_odor[i,:,:],axis = 1) 
        
        cumsum_proba_rest = np.cumsum(proba_density_rest[i,:,:], axis = 1)
        
        #cumsum_proba_odor_file = os.path.join(main_path, "TS_R_timefreq_" + freq + "_all_cond_by_block_trigs",subj_id,sess_index, "plot_cumsum_proba_odor_" + label + ".eps")
        #plot_cormat(cumsum_proba_odor_file,cumsum_proba_odor)
        
        #cumsum_proba_rest_file = os.path.join(main_path, "TS_R_timefreq_" + freq + "_all_cond_by_block_trigs",subj_id,sess_index, "plot_cumsum_proba_rest_" + label + ".eps")
        #plot_cormat(cumsum_proba_rest_file,cumsum_proba_rest)
        
        
        ### diff
        diff_proba_density = cumsum_proba_odor - cumsum_proba_rest
        
        #print diff_proba_density
        
        #diff_proba_density_file = os.path.join(main_path, "TS_R_timefreq_" + freq + "_all_cond_by_block_trigs",subj_id,sess_index, "plot_diff_proba_density_" + label + ".eps")
        #plot_cormat(diff_proba_density_file,diff_proba_density)
        
        list_diff_cumsum_proba.append(diff_proba_density)
        
    np_list_diff_cumsum_proba = np.array(list_diff_cumsum_proba)
    
    print np_list_diff_cumsum_proba.shape
    
    ### by mod
    for i,mod_index in enumerate(np.unique(community_vect)):
            
        plot_mean_mod_diff_probal_density_file = os.path.join(base_dir,subj_id,sess_index, "plot_diff_cumsum_proba_odor-rest_mod_" + freq_modules + "_" + str(mod_index) + ".eps")
        
        print plot_mean_mod_diff_probal_density_file
        
        mean_mod_diff_probal_density = np.mean(np_list_diff_cumsum_proba[community_vect == mod_index,:,:],axis = 0)
        
        
        print mean_mod_diff_probal_density.shape
        
        plot_cormat(plot_mean_mod_diff_probal_density_file,mean_mod_diff_probal_density)
        
def test_sum_proba_density():

    #gather_proba_density_results(subj_id = 'SEMC', freq = '1-100Hz', sess_index = 'R1', trial_cond = 'all')
    
    sum_modular_proba_density(subj_id = 'SEMC', freq = '1-150Hz', freq_modules = 'gamma', sess_index = 'R1', trial_cond = 'all')
    
def test_gather_proba_density():

    gather_proba_density(subj_id = 'VACJ', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
    
    
def test_gather_proba_density_by_chunks():

    #gather_proba_density_by_chunks(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
    #gather_proba_density_only_GM(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
    gather_proba_density_by_description(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
    
def compute_gather_proba_density():

    #for subject_id in subject_ids:
    for subject_id in ['SEMC']:
        for sess_index in ['R1','R2','R3']:
            for trial_cond in ['all']:
                
                gather_proba_density(subj_id = subject_id, freq = '1-150Hz', sess_index = sess_index, trial_cond = trial_cond)
                

def compute_gather_proba_density_by_chunks():

    #for subject_id in subject_ids:
    for subject_id in ['SEMC']:
        for sess_index in ['R1','R2','R3']:
            for trial_cond in ['all']:
                            
                #gather_proba_density_by_chunks(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
                #gather_proba_density_only_GM(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
                gather_proba_density_by_description(subj_id = 'SEMC', freq = '1-150Hz', sess_index = 'R1', trial_cond = 'all')
                
                gather_proba_density_by_chunks(subj_id = subject_id, freq = '1-150Hz', sess_index = sess_index, trial_cond = trial_cond)
                
                
def compute_sum_proba_density():

    for subject_id in subject_ids:
        
        #for freq_modules in freq_band_names:
            
            #for sess_index in ['R1','R2','R3']:
                
                #for trial_cond in ['all']:
                
                    #sum_modular_proba_density(subj_id = subject_id,  freq ='1-150Hz' , freq_modules = freq_modules, sess_index = sess_index, trial_cond = trial_cond)
                 
    
        for sess_index in ['R1','R2','R3']:
            
            for trial_cond in ['all']:
            
                sum_modular_proba_density(subj_id = subject_id,  freq ='1-150Hz' , freq_modules = "beta", sess_index = sess_index, trial_cond = trial_cond)
                
#def compute_sum_modular_ampl():

    #for subject_id in subject_ids:
        
        #for freq in freq_band_names:
            
            #for sess_index in ['R1','R2','R3']:
                
                #for trial_cond in ['all']:
                
                    #sum_modular_ampl(subject_id,  freq, sess_index, trial_cond )
                    
if __name__ == '__main__':
    
    ## run pipeline:
    
        
    #split_correl_analysis_name = correl_analysis_name.split('_')

    #if 'full' in split_correl_analysis_name:
        #main_workflow = create_main_workflow_correct_ampl_by_band()
        
    #elif 'bip' in split_correl_analysis_name:
        
        #if 'odor' in split_correl_analysis_name:
            #main_workflow = create_main_workflow_bip_ampl_by_odor_trigs_by_band()
            
        #elif 'rest' in split_correl_analysis_name:
            #main_workflow = create_main_workflow_bip_ampl_by_rest_trigs_by_band()
    
    #else:
            
        #if 'odor' in split_correl_analysis_name:
            
            #if 'sep' in  split_correl_analysis_name:
                #main_workflow = create_main_workflow_ampl_by_sep_odor_trigs_by_band()
                
            #else:
                #main_workflow = create_main_workflow_ampl_by_odor_trigs_by_band()
                
        #elif 'rest' in split_correl_analysis_name:
            
            #if 'sep' in  split_correl_analysis_name:
                #main_workflow = create_main_workflow_ampl_by_sep_rest_trigs_by_band()
                
            #else:
                #main_workflow = create_main_workflow_ampl_by_rest_trigs_by_band()
                
    ######## run
    #main_workflow.config['execution'] = {'remove_unnecessary_outputs':'false'}
    #main_workflow.run(plugin='MultiProc', plugin_args={'n_procs' : 8})    

    ########################################## gather results #################################
    
    #test_gather_proba_density()
    #compute_gather_proba_density()
    
    # by chunks
    #test_gather_proba_density_by_chunks()
    compute_gather_proba_density_by_chunks()
    
    ########################################## sum and plot resulting signals by modules ##############
    
    #test_sum_proba_density()
    #compute_sum_proba_density()
    
    
    
    
    
    