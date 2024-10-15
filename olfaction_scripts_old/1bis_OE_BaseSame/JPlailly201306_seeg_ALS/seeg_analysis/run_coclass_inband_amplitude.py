# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 

from connection import *

from params import main_path,joblib_cachedir
from params import subject_ids,exp,trial_R_conds,trial_E_conds

from params import conf_interval_prob,Z_thr
from params import radatools_optim,radatools_path

from params import freq_band_names

from params import coclass_analysis_name,split_coclass_analysis_name

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

from dmgraphanalysis_nodes.nodes.coclass import PrepareCoclass,PlotCoclass,PlotIGraphCoclass,DiffMatrices,PlotIGraphConjCoclass
from dmgraphanalysis_nodes.nodes.modularity import ComputeIntNetList, PrepRada, CommRada, PlotIGraphModules, ComputeNodeRoles, NetPropRada

from dmgraphanalysis_nodes.utils import get_first

#def create_infosource_R_correct_cond_ampl_by_band():

    #infosource = pe.Node(interface=IdentityInterface(fields=['subject_id', 'sess_index','trial_cond','freq_band']),name="infosource")
    #infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1','R2','R3']),('trial_cond',trial_R_conds),('freq_band',freq_band_names)]

    #### test 
    ##infosource.iterables = [('subject_id', ['SEMC']),('sess_index',['R1','R2','R3']),('trial_cond',['all']),('freq_band',['beta'])]
    ##infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1']),('trial_cond',['all']),('freq_band',['beta'])]
    
    #return infosource
    
def create_infosource_R_correct_cond_trial_ampl_by_band():

    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id', 'sess_index','cond','trial_cond','freq_band']),name="infosource")
    infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1','R2','R3']),('cond',['odor','rest']),('trial_cond',['all']),('freq_band',freq_band_names)]

    ### test 
    #infosource.iterables = [('subject_id', ['SEMC']),('sess_index',['R1','R2','R3']),('trial_cond',['all']),('freq_band',['beta'])]
    #infosource.iterables = [('subject_id', subject_ids),('sess_index',['R1']),('trial_cond',['all']),('freq_band',['beta'])]
    
    return infosource
    
def create_datasource_ampl_by_trigs_by_band():
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','sess_index','cond','trial_cond','freq_band'],outfields=['mod_files','node_corres_files','coords_files']),name = 'datasource')
    
    datasource.inputs.base_directory = main_path
    
    datasource.inputs.template = '%s%s%s/%s%s%s%s%s%s%s%s/%s/%s/%s*/%s'
    datasource.inputs.template_args = dict(
    
    mod_files = [["correl_ampl_by_sep_",'cond',"_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"),"_freq_band_",'freq_band',"_sess_index_",'sess_index',"_subject_id_",'subject_id',"_trial_cond_",'trial_cond',"community_rada","mapflow","_community_rada","Z_List.lol"]],
    node_corres_files = [["correl_ampl_by_sep_",'cond',"_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"),"_freq_band_",'freq_band',"_sess_index_",'sess_index',"_subject_id_",'subject_id',"_trial_cond_",'trial_cond',"prep_rada","mapflow","_prep_rada","Z_List.net"]],
    coords_files = [["correl_ampl_by_sep_",'cond',"_trigs_by_band-Z_thr_" + str(Z_thr).replace(".","_"),"_freq_band_",'freq_band',"_sess_index_",'sess_index',"_subject_id_",'subject_id',"_trial_cond_",'trial_cond',"compute_net_List","mapflow","_compute_net_List","coords.txt"]],
        )
    datasource.inputs.sort_filelist = True
    return datasource
    
    
def create_main_workflow_coclass_by_sep_trigs_by_band():

    main_workflow = pe.Workflow(name=coclass_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    if exp == 'R':
        infosource = create_infosource_R_correct_cond_trial_ampl_by_band()
        
    elif exp == 'E':
        
        print "Not implemented yet"
        
        sys.exit()
        
    
    ## Data source
    datasource = create_datasource_ampl_by_trigs_by_band()
    
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')
    main_workflow.connect(infosource, 'cond', datasource, 'cond')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'freq_band', datasource, 'freq_band')
    
    ##################################### compute sum coclass and group based coclass matrices  #####################################################
    
    #### prepare_nbs_stats_rada
    
    prepare_coclass = pe.Node(interface = PrepareCoclass(),name='prepare_coclass')
    
    
    main_workflow.connect(datasource, 'mod_files',prepare_coclass,'mod_files')
    main_workflow.connect(datasource, 'node_corres_files',prepare_coclass,'node_corres_files')
    main_workflow.connect(datasource, 'coords_files',prepare_coclass,'coords_files')
    
    main_workflow.connect(datasource, ('coords_files',get_first),prepare_coclass,'gm_mask_coords_file')
    
    ######################################################## coclassification matrices ############################################################
    
    ######### norm coclass
    plot_norm_coclass= pe.Node(interface = PlotCoclass(),name='plot_norm_coclass')
    
    #plot_norm_coclass.inputs.labels_file = ROI_coords_labels_file
    plot_norm_coclass.inputs.list_value_range = [0,100]
    
    main_workflow.connect(prepare_coclass, 'norm_coclass_matrix_file',plot_norm_coclass,'coclass_matrix_file')
    
    
    plot_igraph_norm_coclass= pe.Node(interface = PlotIGraphCoclass(),name='plot_igraph_norm_coclass')
    #plot_igraph_norm_coclass.inputs.gm_mask_coords_file = ROI_coords_MNI_coords_file
    plot_igraph_norm_coclass.inputs.threshold = 50
    #plot_igraph_norm_coclass.inputs.labels_file = ROI_coords_labels_file
    
    main_workflow.connect(prepare_coclass, 'norm_coclass_matrix_file',plot_igraph_norm_coclass,'coclass_matrix_file')
    main_workflow.connect(datasource, ('coords_files',get_first),plot_igraph_norm_coclass,'gm_mask_coords_file')
    
    ################################################################ modular decomposition on norm_coclass ######################################################################################
    
    ########### serie 
    
    if 'rada' in split_coclass_analysis_name:
        
        ### compute Z_list from coclass matrix
        compute_list_norm_coclass = pe.Node(interface = ComputeIntNetList(),name='compute_list_norm_coclass')
        compute_list_norm_coclass.inputs.threshold = 50
        
        main_workflow.connect(prepare_coclass,'norm_coclass_matrix_file',compute_list_norm_coclass, 'int_mat_file')
        
        
        #################################################### radatools ################################################################

        ### prepare net_list for radatools processing  
        prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
        prep_rada.inputs.radatools_path = radatools_path
        
        main_workflow.connect(compute_list_norm_coclass, 'net_List_file', prep_rada, 'net_List_file')
        
        ### compute community with radatools
        community_rada = pe.Node(interface = CommRada(), name='community_rada')
        community_rada.inputs.optim_seq = radatools_optim
        community_rada.inputs.radatools_path = radatools_path
        
        main_workflow.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
        
        ### node roles
        node_roles = pe.Node(interface = ComputeNodeRoles(role_type = "4roles"), name='node_roles')
        
        main_workflow.connect( prep_rada, 'Pajek_net_file',node_roles,'Pajek_net_file')
        main_workflow.connect( community_rada, 'rada_lol_file',node_roles,'rada_lol_file')
        

            
        #### plot_igraph_modules_rada_norm_coclass
        
        plot_igraph_modules_rada_norm_coclass = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada_norm_coclass')
        
        #plot_igraph_modules_rada_norm_coclass.inputs.labels_file = ROI_coords_labels_file
        #plot_igraph_modules_rada_norm_coclass.inputs.coords_file = ROI_coords_MNI_coords_file
        
        main_workflow.connect(datasource, ('coords_files',get_first),plot_igraph_modules_rada_norm_coclass,'coords_file')
    
        main_workflow.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada_norm_coclass,'Pajek_net_file')
        main_workflow.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada_norm_coclass,'rada_lol_file')
        
        main_workflow.connect(node_roles, 'node_roles_file',plot_igraph_modules_rada_norm_coclass,'node_roles_file')
        
        
        
        ############ compute networ properties with rada
        
        net_prop = pe.Node(interface = NetPropRada(optim_seq = "A"), name = 'net_prop')
        net_prop.inputs.radatools_path = radatools_path
        
        main_workflow.connect(prep_rada, 'Pajek_net_file',net_prop,'Pajek_net_file')
        
    return main_workflow

if __name__ == '__main__':
    
    ### run pipeline:
    main_workflow = create_main_workflow_coclass_by_sep_trigs_by_band()
    
    
    ###### run
    main_workflow.config['execution'] = {'remove_unnecessary_outputs':'false'}
    main_workflow.run(plugin='MultiProc', plugin_args={'n_procs' : 8})    

    
    ########################################## gather results #################################
    
    
    