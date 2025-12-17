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

#t_win_start,t_win_stop,baseline_t_win_start,baseline_t_win_stop,f_start,f_stop,envelop_method,baseline_mode

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

def create_infosource_R_correct_ampl_by_band():

    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id', 'sess_index','freq_band']),name="infosource")
    infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1','R2','R3']),('freq_band',freq_band_names)]

    ### test 
    #infosource.iterables = [('subject_id', ['CHAF']),('sess_index',['R1']),('freq_band',['beta'])]
    
    return infosource
    
def create_infosource_R_correct_cond_ampl_by_band():

    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id', 'sess_index','trial_cond','freq_band']),name="infosource")
    infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1','R2','R3']),('trial_cond',trial_R_conds),('freq_band',freq_band_names)]
    
    #infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1','R2','R3']),('trial_cond',['all']),('freq_band',freq_band_names)]

    ### test
    #infosource.iterables = [('subject_id', ['SEMC']),('sess_index',['R1']),('trial_cond',['all']),('freq_band',['gamma'])]
    
    #infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1']),('trial_cond',['all']),('freq_band',['beta'])]
    
    return infosource
    
########################################################################## all trig ampl #################################################################################

def create_datasource_correct_ampl_by_band():
        
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','freq_band'],outfields=['ts_trig_ampl_file','channel_names_file','channel_coords_file']),name = 'datasource')
    datasource.inputs.base_directory = main_path
    datasource.inputs.template = '%s%s%s/%s/%s/%s'
    datasource.inputs.template_args = dict(
        ts_trig_ampl_file=[["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id','sess_index',"np_correct_ampl_sigs.npy"]],
        channel_names_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_names.txt"]],
        channel_coords_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_coords.txt"]]
        )
    datasource.inputs.sort_filelist = True
    return datasource
    
def create_main_workflow_correct_ampl_by_band():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    infosource = create_infosource_R_correct_ampl_by_band()
    
    ## Data source
    datasource = create_datasource_correct_ampl_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ### correl_mat
    compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    main_workflow.connect(datasource, 'ts_trig_ampl_file', compute_conf_cor_mat, 'ts_file')
    main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    #### net_list
    compute_net_List = pe.Node(interface = ComputeNetList(threshold = Z_thr),name='compute_net_List')
    
    main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    ##### net_list
    #compute_net_List = pe.Node(interface = ComputeNetList(),name='compute_net_List')
    
    #main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    #main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    #################################################### radatools ################################################################

    ### prepare net_list for radatools processing  
    prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    #Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    ### compute community with radatools
    community_rada = pe.Node(interface = CommRada(), name='community_rada')
    community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    ### plot radatools
    
    #### plot_igraph_modules_rada
    
    plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
    
    main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    return main_workflow
    
####################################################################### by trigs ############################################################################################################

def create_datasource_ampl_by_odor_trigs_by_band():
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','trial_cond','freq_band'],outfields=['ts_ampl_by_odor_trig_file','channel_names_file','channel_coords_file']),name = 'datasource')
    datasource.inputs.base_directory = main_path
    datasource.inputs.template = '%s%s%s/%s/%s/%s%s%s'
    datasource.inputs.template_args = dict(
        ts_ampl_by_odor_trig_file=[["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id','sess_index',"correct_ts_",'trial_cond',"_ampl_by_odor_trigs.npy"]],
        channel_names_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_names.txt","",""]],
        channel_coords_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_coords.txt","",""]]
        )
    datasource.inputs.sort_filelist = True
    return datasource
    
def create_main_workflow_ampl_by_odor_trigs_by_band():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    if exp == 'R':
        infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    elif exp == 'E':
        
        print "Not implemented yet"
        
        sys.exit()
        
    
    ## Data source
    datasource = create_datasource_ampl_by_odor_trigs_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ### concat cond ts
    concat_ts = pe.Node(interface = ConcatTS(),name = 'concat_ts')
    
    main_workflow.connect(datasource, 'ts_ampl_by_odor_trig_file', concat_ts, 'all_ts_file')
    
    
    ### correl_mat
    compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    main_workflow.connect(concat_ts, 'concatenated_ts_file', compute_conf_cor_mat, 'ts_file')
    main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    #### net_list
    compute_net_List = pe.Node(interface = ComputeNetList(threshold = Z_thr),name='compute_net_List')
    
    main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    #################################################### radatools ################################################################

    ### prepare net_list for radatools processing  
    prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    #Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    ### compute community with radatools
    community_rada = pe.Node(interface = CommRada(), name='community_rada')
    community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    ### plot radatools
    
    #### plot_igraph_modules_rada
    
    plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
    
    main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    return main_workflow


def create_main_workflow_ampl_by_sep_odor_trigs_by_band():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    
    main_workflow.base_dir = main_path
    
    ## Info source
    if exp == 'R':
        infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    elif exp == 'E':
        
        print "Not implemented yet"
        
        sys.exit()
        
    ## Data source
    datasource = create_datasource_ampl_by_odor_trigs_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ### concat cond ts
    sep_ts = pe.Node(interface = SeparateTS(),name = 'sep_ts')
    
    main_workflow.connect(datasource, 'ts_ampl_by_odor_trig_file', sep_ts, 'all_ts_file')
    
    
    ### correl_mat
    compute_conf_cor_mat = pe.MapNode(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    main_workflow.connect(sep_ts, 'separated_ts_files', compute_conf_cor_mat, 'ts_file')
    main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    #### net_list
    compute_net_List = pe.MapNode(interface = ComputeNetList(threshold = Z_thr),iterfield = ['Z_cor_mat_file'],name='compute_net_List')
    
    main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    #################################################### radatools ################################################################

    ### prepare net_list for radatools processing  
    prep_rada = pe.MapNode(interface = PrepRada(),iterfield = ['net_List_file'],name='prep_rada')
    #Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    ### compute community with radatools
    community_rada = pe.MapNode(interface = CommRada(),iterfield = ['Pajek_net_file'], name='community_rada')
    community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    ### plot radatools
    
    #### plot_igraph_modules_rada
    
    plot_igraph_modules_rada = pe.MapNode(interface = PlotIGraphModules(),iterfield = ['Pajek_net_file','rada_lol_file'],name='plot_igraph_modules_rada')
    
    main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    return main_workflow

#### rest
def create_datasource_ampl_by_rest_trigs_by_band():
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','trial_cond','freq_band'],outfields=['ts_ampl_by_rest_trig_file','channel_names_file','channel_coords_file']),name = 'datasource')
    datasource.inputs.base_directory = main_path
    datasource.inputs.template = '%s%s%s/%s/%s/%s%s%s'
    datasource.inputs.template_args = dict(
        ts_ampl_by_rest_trig_file=[["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id','sess_index',"correct_ts_",'trial_cond',"_ampl_by_rest_trigs.npy"]],
        channel_names_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_names.txt","",""]],
        channel_coords_file = [["TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_channel_coords.txt","",""]]
        )
    datasource.inputs.sort_filelist = True
    return datasource
    
def create_main_workflow_ampl_by_rest_trigs_by_band():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    if exp == 'R':
        infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    elif exp == 'E':
        
        print "Not implemented yet"
        
        sys.exit()
        
    
    ## Data source
    datasource = create_datasource_ampl_by_rest_trigs_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ### concat cond ts
    concat_ts = pe.Node(interface = ConcatTS(),name = 'concat_ts')
    
    main_workflow.connect(datasource, 'ts_ampl_by_rest_trig_file', concat_ts, 'all_ts_file')
    
    
    ### correl_mat
    compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    main_workflow.connect(concat_ts, 'concatenated_ts_file', compute_conf_cor_mat, 'ts_file')
    main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    #### net_list
    compute_net_List = pe.Node(interface = ComputeNetList(threshold = Z_thr),name='compute_net_List')
    
    main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    #################################################### radatools ################################################################

    ### prepare net_list for radatools processing  
    prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    #Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    ### compute community with radatools
    community_rada = pe.Node(interface = CommRada(), name='community_rada')
    community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    ### plot radatools
    
    #### plot_igraph_modules_rada
    
    plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
    
    main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    return main_workflow
    
def create_main_workflow_ampl_by_sep_rest_trigs_by_band():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    if exp == 'R':
        infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    elif exp == 'E':
        
        print "Not implemented yet"
        
        sys.exit()
        
    
    ## Data source
    datasource = create_datasource_ampl_by_rest_trigs_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ### concat cond ts
    sep_ts = pe.Node(interface = SeparateTS(),name = 'sep_ts')
    
    main_workflow.connect(datasource, 'ts_ampl_by_rest_trig_file', sep_ts, 'all_ts_file')
    
    
    ### correl_mat
    compute_conf_cor_mat = pe.MapNode(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    main_workflow.connect(sep_ts, 'separated_ts_files', compute_conf_cor_mat, 'ts_file')
    main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    #### net_list
    compute_net_List = pe.MapNode(interface = ComputeNetList(threshold = Z_thr),iterfield = ['Z_cor_mat_file'],name='compute_net_List')
    
    main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    #################################################### radatools ################################################################

    ### prepare net_list for radatools processing  
    prep_rada = pe.MapNode(interface = PrepRada(),iterfield = ['net_List_file'],name='prep_rada')
    #Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    ### compute community with radatools
    community_rada = pe.MapNode(interface = CommRada(),iterfield = ['Pajek_net_file'], name='community_rada')
    community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    ### plot radatools
    
    #### plot_igraph_modules_rada
    
    plot_igraph_modules_rada = pe.MapNode(interface = PlotIGraphModules(),iterfield = ['Pajek_net_file','rada_lol_file'],name='plot_igraph_modules_rada')
    
    main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    return main_workflow

    
    
    
    
def gather_cormat_results_ampl_by_odor_and_rest():
    
    from dmgraphanalysis_nodes.utils_plot import plot_ranged_cormat,plot_int_mat
    
    from dmgraphanalysis_nodes.utils_stats import compute_pairwise_ttest_fdr
    import glob
    
    
    
    from plot_electrodes import plot_electrodes_graphs
    
    #for subj_id in subject_ids:
    for subj_id in ['SEMC']:
      
        print subj_id
        
        for freq in freq_band_names:
            
            base_dir = os.path.join(main_path,"TS_" + exp + "_" + envelop_method + "_"+ freq +"_all_cond_by_block_trigs")

            print freq
            
            labels_file = os.path.join(base_dir,subj_id ,"correct_channel_names.txt")
        
            labels = [line.strip() for line in open(labels_file)]
            print len(labels)
            
            coords_file = os.path.join(base_dir,subj_id ,"correct_channel_coords.txt")
            
            coords = np.array(np.loadtxt(coords_file),dtype = float)
            
            for trial in trial_R_conds:
                
                print coords
            
                odor_cormat_files = glob.glob(os.path.join(main_path,"correl_ampl_by_sep_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"),"_freq_band_" + freq + "_sess_index_*_subject_id_" + subj_id + "_trial_cond_" + trial ,"compute_conf_cor_mat","mapflow","_compute_conf_cor_mat*","Z_cor_mat_correct_ts_" + trial + "_ampl_by_odor_trigs_trig_*.npy"))
                
                print odor_cormat_files
                
                rest_cormat_files = glob.glob(os.path.join(main_path,"correl_ampl_by_sep_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"),"_freq_band_" + freq + "_sess_index_*_subject_id_" + subj_id + "_trial_cond_" + trial,"compute_conf_cor_mat","mapflow","_compute_conf_cor_mat*","Z_cor_mat_correct_ts_" + trial + "_ampl_by_rest_trigs_trig_*.npy"))
                
                print rest_cormat_files
                
                print len(odor_cormat_files),len(rest_cormat_files)
                
                if len(odor_cormat_files) == len(rest_cormat_files) and 5 <= len(odor_cormat_files):
                    
                    ### odor
                    group_odor_cormat = []
                    
                    for i,odor_cormat_file in enumerate(odor_cormat_files):
                        
                        odor_cormat = np.load(odor_cormat_file)
                        
                        group_odor_cormat.append(odor_cormat + np.transpose(odor_cormat))
                    
                    group_odor_cormat = np.array(group_odor_cormat,dtype = 'float')
                    
                    new_group_odor_cormat = np.swapaxes(group_odor_cormat,0,2)
                    
                    print new_group_odor_cormat.shape
                    
                    ### rest
                    group_rest_cormat = []
                    
                    for i,rest_cormat_file in enumerate(rest_cormat_files):
                        
                        rest_cormat = np.load(rest_cormat_file)
                        
                        group_rest_cormat.append(rest_cormat + np.transpose(rest_cormat))
                    
                    group_rest_cormat = np.array(group_rest_cormat,dtype = 'float')
                    
                    new_group_rest_cormat = np.swapaxes(group_rest_cormat,0,2)
                    
                    print new_group_rest_cormat.shape
                    
                    #### differences
                    
                    diff_odor_rest_mat = new_group_odor_cormat - new_group_rest_cormat
                    
                    print diff_odor_rest_mat.shape
                    
                    diff_odor_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_odor_rest_" + freq + "_" + trial + "_subject_id_" + subj_id + ".npy")
                    
                    np.save(diff_odor_rest_mat_file,diff_odor_rest_mat)
                    
                    mean_diff_odor_rest_mat = np.mean(diff_odor_rest_mat,axis = 2)
                    
                    print mean_diff_odor_rest_mat.shape
                    
                    plot_mean_diff_odor_rest_mat_file = os.path.join(main_path,correl_analysis_name,"mean_diff_odor_rest_" + freq + "_"+ trial + "_subject_id_" + subj_id + ".eps")
                    
                    plot_ranged_cormat(plot_mean_diff_odor_rest_mat_file,mean_diff_odor_rest_mat,list_labels = labels,fix_full_range = [-1.0,1.0])
                                
                    #### signif differences
                    diff_signif_odor_rest_mat = compute_pairwise_ttest_fdr(new_group_odor_cormat,new_group_rest_cormat,t_test_thresh_fdr = 0.05,paired = True)
                    
                    diff_signif_odor_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_odor_rest_" + freq + "_"+ trial + "_subject_id_" + subj_id + ".npy")
                    
                    np.save(diff_signif_odor_rest_mat_file,diff_signif_odor_rest_mat)
                    
                    plot_diff_signif_odor_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_odor_rest_" + freq + "_"+ trial + "_subject_id_" + subj_id + ".eps")
                    
                    plot_int_mat(plot_diff_signif_odor_rest_mat_file,diff_signif_odor_rest_mat,list_labels = labels,fix_full_range = [-4.0,4.0])
                
                ## plotting with pyplotbrain
                plot_electrodes_graphs(labels,coords,diff_signif_odor_rest_mat,export_path = os.path.join(main_path,correl_analysis_name),pref = "diff_signif_odor_rest_"+ freq + "_" + trial + "_subject_id_" + subj_id )
                
        ########## testing correlation differences (odor-rest) between conditions
        
        #print "testing correlation differences (odor-rest) between conditions"
        
        #diff_hit_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_odor_rest_hit_subject_id_" + subj_id + ".npy")
        #diff_fa_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_odor_rest_fa_subject_id_" + subj_id + ".npy")
        #diff_cr_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_odor_rest_cr_subject_id_" + subj_id + ".npy")
        
        #if os.path.exists(diff_hit_rest_mat_file) and os.path.exists(diff_fa_rest_mat_file):
        
            #print "Ok for hit-rest vs fa-rest"
            
            #diff_hit_rest_mat = np.load(diff_hit_rest_mat_file)
            
            #print diff_hit_rest_mat.shape
            
            #diff_fa_rest_mat = np.load(diff_fa_rest_mat_file)
            
            #print diff_fa_rest_mat.shape
            
            
            ##### signif differences
            #diff_signif_hit_rest_fa_rest_mat = compute_pairwise_ttest_fdr(diff_hit_rest_mat,diff_fa_rest_mat,t_test_thresh_fdr = 0.05,paired = True)
            
            #diff_signif_hit_rest_fa_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_hit-rest_fa-rest_subject_id_" + subj_id + ".npy")
            
            #np.save(diff_signif_hit_rest_fa_rest_mat_file,diff_signif_hit_rest_fa_rest_mat)
            
            #plot_diff_signif_hit_rest_fa_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_hit-rest_fa-rest_subject_id_" + subj_id + ".eps")
            
            #plot_int_mat(plot_diff_signif_hit_rest_fa_rest_mat_file,diff_signif_hit_rest_fa_rest_mat,list_labels = labels,fix_full_range = [-4.0,4.0])
        
            ### plotting with pyplotbrain
            #plot_electrodes_graphs(labels,coords,diff_signif_hit_rest_fa_rest_mat,export_path = os.path.join(main_path,correl_analysis_name),pref = "diff_signif_hit-rest_fa-rest_subject_id_" + subj_id )
            
        #if os.path.exists(diff_hit_rest_mat_file) and os.path.exists(diff_cr_rest_mat_file):
        
            #print "Ok for hit-rest vs cr-rest "
            
            #diff_hit_rest_mat = np.load(diff_hit_rest_mat_file)
            
            #print diff_hit_rest_mat.shape
            
            #diff_cr_rest_mat = np.load(diff_cr_rest_mat_file)
            
            #print diff_cr_rest_mat.shape
            
            
            ##### signif differences
            #diff_signif_hit_rest_cr_rest_mat = compute_pairwise_ttest_fdr(diff_hit_rest_mat,diff_cr_rest_mat,t_test_thresh_fdr = 0.05,paired = False)
            
            #diff_signif_hit_rest_cr_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_hit-rest_cr-rest_subject_id_" + subj_id + ".npy")
            
            #np.save(diff_signif_hit_rest_cr_rest_mat_file,diff_signif_hit_rest_cr_rest_mat)
            
            #plot_diff_signif_hit_rest_cr_rest_mat_file = os.path.join(main_path,correl_analysis_name,"diff_signif_hit-rest_cr-rest_subject_id_" + subj_id + ".eps")
            
            #plot_int_mat(plot_diff_signif_hit_rest_cr_rest_mat_file,diff_signif_hit_rest_cr_rest_mat,list_labels = labels,fix_full_range = [-4.0,4.0])
        
            #### plotting with pyplotbrain
            #plot_electrodes_graphs(labels,coords,diff_signif_hit_rest_cr_rest_mat,export_path = os.path.join(main_path,correl_analysis_name),pref = "diff_signif_hit-rest_cr-rest_subject_id_" + subj_id )
            
##################################################################### by trigs bipolar ##########################################################


#def create_datasource_bip_ampl_by_odor_trigs_by_band():
    
    #datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','trial_cond','freq_band'],outfields=['ts_ampl_by_odor_trig_file','channel_names_file','channel_coords_file']),name = 'datasource')
    #datasource.inputs.base_directory = main_path
    #datasource.inputs.template = '%s%s%s/%s/%s/%s%s%s'
    #datasource.inputs.template_args = dict(
        #ts_ampl_by_odor_trig_file=[["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id','sess_index',"correct_ts_",'trial_cond',"_ampl_by_odor_trigs.npy"]],
        #channel_names_file = [["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_bip_channel_names.txt","",""]],
        #channel_coords_file = [["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_bip_channel_coords.txt","",""]]
        #)
    #datasource.inputs.sort_filelist = True
    #return datasource
    
#def create_main_workflow_bip_ampl_by_odor_trigs_by_band():

    #main_workflow = pe.Workflow(name="correl_bip_ampl_by_odor_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"))
    #main_workflow.base_dir = main_path
    
    ### Info source
    #if exp == 'R':
        #infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    #elif exp == 'E':
        
        #print "Not implemented yet"
        
        #sys.exit()
        
    
    ### Data source
    #datasource = create_datasource_bip_ampl_by_odor_trigs_by_band()
    
    #main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    #main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    #main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    #main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    #### concat cond ts
    #concat_ts = pe.Node(interface = ConcatTS(),name = 'concat_ts')
    
    #main_workflow.connect(datasource, 'ts_ampl_by_odor_trig_file', concat_ts, 'all_ts_file')
    
    
    #### correl_mat
    #compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    #compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    #main_workflow.connect(concat_ts, 'concatenated_ts_file', compute_conf_cor_mat, 'ts_file')
    #main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    ##### net_list
    #compute_net_List = pe.Node(interface = ComputeNetList(threshold = Z_thr),name='compute_net_List')
    
    #main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    #main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    ##################################################### radatools ################################################################

    #### prepare net_list for radatools processing  
    #prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    ##Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    #prep_rada.inputs.radatools_path = radatools_path
    
    #main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    #### compute community with radatools
    #community_rada = pe.Node(interface = CommRada(), name='community_rada')
    #community_rada.inputs.optim_seq = radatools_optim
    #community_rada.inputs.radatools_path = radatools_path
    
    #main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    #### plot radatools
    
    ##### plot_igraph_modules_rada
    
    #plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
    
    #main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    #main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    #main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    #main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    #return main_workflow

#def create_datasource_bip_ampl_by_rest_trigs_by_band():
    
    #datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','trial_cond','freq_band'],outfields=['ts_ampl_by_rest_trig_file','channel_names_file','channel_coords_file']),name = 'datasource')
    #datasource.inputs.base_directory = main_path
    #datasource.inputs.template = '%s%s%s/%s/%s/%s%s%s'
    #datasource.inputs.template_args = dict(
        #ts_ampl_by_rest_trig_file=[["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id','sess_index',"correct_ts_",'trial_cond',"_ampl_by_rest_trigs.npy"]],
        #channel_names_file = [["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_bip_channel_names.txt","",""]],
        #channel_coords_file = [["bip_TS_R_tfr_",'freq_band',"_all_cond_by_block_trigs",'subject_id',"","correct_bip_channel_coords.txt","",""]]
        #)
    #datasource.inputs.sort_filelist = True
    #return datasource
    
#def create_main_workflow_bip_ampl_by_rest_trigs_by_band():

    #main_workflow = pe.Workflow(name="correl_bip_ampl_by_rest_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"))
    #main_workflow.base_dir = main_path
    
    ### Info source
    #if exp == 'R':
        #infosource = create_infosource_R_correct_cond_ampl_by_band()
        
    #elif exp == 'E':
        
        #print "Not implemented yet"
        
        #sys.exit()
        
    
    ### Data source
    #datasource = create_datasource_bip_ampl_by_rest_trigs_by_band()
    
    #main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    #main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    #main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    #main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    #### concat cond ts
    #concat_ts = pe.Node(interface = ConcatTS(),name = 'concat_ts')
    
    #main_workflow.connect(datasource, 'ts_ampl_by_rest_trig_file', concat_ts, 'all_ts_file')
    
    
    #### correl_mat
    #compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),iterfield = ['ts_file'],name='compute_conf_cor_mat')
    
    #compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    #main_workflow.connect(concat_ts, 'concatenated_ts_file', compute_conf_cor_mat, 'ts_file')
    #main_workflow.connect(datasource, 'channel_names_file', compute_conf_cor_mat, 'labels_file')
    
    ##### net_list
    #compute_net_List = pe.Node(interface = ComputeNetList(threshold = Z_thr),name='compute_net_List')
    
    #main_workflow.connect(compute_conf_cor_mat, 'Z_cor_mat_file',compute_net_List, 'Z_cor_mat_file')
    #main_workflow.connect(datasource, 'channel_coords_file',compute_net_List, 'coords_file')
    
    
    ##################################################### radatools ################################################################

    #### prepare net_list for radatools processing  
    #prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    ##Function(input_names=['List_net_file','radatools_prep_path'],output_names = ['Pajek_net_file'],function = prep_radatools),name='prep_rada')
    #prep_rada.inputs.radatools_path = radatools_path
    
    #main_workflow.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
    #### compute community with radatools
    #community_rada = pe.Node(interface = CommRada(), name='community_rada')
    #community_rada.inputs.optim_seq = radatools_optim
    #community_rada.inputs.radatools_path = radatools_path
    
    #main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    
    #### plot radatools
    
    ##### plot_igraph_modules_rada
    
    #plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
    
    #main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
    #main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
    
    #main_workflow.connect(datasource, 'channel_coords_file',plot_igraph_modules_rada,'coords_file')
    #main_workflow.connect(datasource, 'channel_names_file',plot_igraph_modules_rada,'labels_file')
    
    #return main_workflow
    
def sum_modular_ampl(subj_id = 'SEMC', freq = 'beta', sess_index = 'R1', trial_cond = 'all'):

    from dmgraphanalysis_nodes.utils_net import read_lol_file

    from dmgraphanalysis_nodes.utils_plot import plot_signals
    
    from dmgraphanalysis_nodes.utils_igraph import igraph_colors
    
    from plot_electrodes import hex_to_rgb
    
    
    path = os.path.join(main_path, correl_analysis_name, "_freq_band_" + freq + "_sess_index_" + sess_index + "_subject_id_" + subj_id + "_trial_cond_" + trial_cond)
    
    lol_file = os.path.join(path,"community_rada","Z_List.lol")
    
    community_vect = read_lol_file(lol_file)
    
    print community_vect
    
    ######## concat sig
    #concatenated_ts_file = os.path.join(path,"concat_ts","concatenated_ts.npy")
    
    #concatenated_ts = np.load(concatenated_ts_file)
    
    #print concatenated_ts.shape
    
    #sum_sig = np.array([np.sum(concatenated_ts[mod_index == community_vect], axis = 0) for mod_index in np.unique(community_vect)])
    
    #print sum_sig.shape
    
    #plot_sig_file = os.path.join(main_path, correl_analysis_name, "plot_concat_sig_by_mod.eps")
    
    #plot_signals(plot_sig_file,sum_sig)
    
    ####### orig sig
    orig_ts_by_trig_file = os.path.join(main_path, "TS_R_tfr_" + freq + "_all_cond_by_block_trigs",subj_id,sess_index,"correct_ts_" + trial_cond + "_ampl_by_odor_trigs.npy")
    
    orig_ts_by_trig = np.load(orig_ts_by_trig_file)
    
    print orig_ts_by_trig.shape
    
    colors_by_mod = [igraph_colors[i] for i in np.unique(community_vect)]
    
    print colors_by_mod
    
    colors_plots_by_mod = np.array(colors_by_mod,dtype = 'str')[community_vect]
    
    print colors_plots_by_mod
        
    #### by trig
    
    for index_trig in range(orig_ts_by_trig.shape[0]):
    
        print index_trig
        
        trig_ts = orig_ts_by_trig[index_trig,:,:]
        
        mean_sig = np.array([np.mean(trig_ts[mod_index == community_vect], axis = 0) for mod_index in np.unique(community_vect)])
        
        print mean_sig.shape
        
        plot_mean_sig_file = os.path.join(path, "plot_mean_sig_odor_trig_" + str(index_trig) + ".eps")
        
        plot_signals(plot_mean_sig_file,mean_sig,colors_by_mod,ylim = [0,10])
        
        #### 
        
        #plot_all_sig_file = os.path.join(path, "plot_all_sig_odor_trig_" + str(index_trig) + ".eps")
        
        #print plot_all_sig_file
        
        #plot_signals(plot_all_sig_file,trig_ts,colors_plots_by_mod,ylim = [0,10])
            
        for i,mod_index in enumerate(np.unique(community_vect)):
                
            plot_all_sig_file = os.path.join(path, "plot_sig_mod_" + str(mod_index) + "_odor_trig_" + str(index_trig) + ".eps")
            
            print plot_all_sig_file
            
            plot_signals(plot_all_sig_file,trig_ts[community_vect == mod_index],[colors_by_mod[i]],ylim = [0,10])
                
    ##### by mod
    
    #for i,mod_index in enumerate(np.unique(community_vect)):
    
        #print mod_index
        
        #mean_sig = np.array([np.mean(orig_ts_by_trig[index_trig,:,:][mod_index == community_vect], axis = 0) for index_trig in range(orig_ts_by_trig.shape[0])])
        
        #print mean_sig.shape
        
        #plot_sig_file = os.path.join(path, "plot_sig_odor_mod_" + str(mod_index) + ".eps")
        
        ##print igraph_colors[mod_index]
                
        #plot_signals(plot_sig_file,mean_sig, colors = [colors_by_mod[i]])
        
        
def test_sum_modular_ampl():

    sum_modular_ampl(subj_id = 'SEMC', freq = 'gamma', sess_index = 'R1', trial_cond = 'all')
    
def compute_sum_modular_ampl():

    for subject_id in subject_ids:
        
        for freq in freq_band_names:
            
            for sess_index in ['R1','R2','R3']:
                
                for trial_cond in ['all']:
                
                    sum_modular_ampl(subject_id,  freq, sess_index, trial_cond )
                    
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
    
    #gather_cormat_results_ampl_by_odor_and_rest()
    
    
    ########################################## sum and plot resulting signals by modules ##############
    
    test_sum_modular_ampl()
    
    #compute_sum_modular_ampl()
    
    
    
    
    