# -*- coding: utf-8 -*-

import sys, os

from define_variables import main_path, data_path
from define_variables import subject_ids

from define_variables import freq_band_names,freq_bands
from define_variables import con_method
#,epoch_window_length
from define_variables import correl_analysis_name

from define_variables import exp

#from define_variables import column_label

from define_variables import sfreq

#from define_variables import con_den

#from define_variables import MEG_elec_coords_file ,MEG_elec_names_file
#,con_thr
    
#from define_variables import exp,envelop_method

#t_win_start,t_win_stop,baseline_t_win_start,baseline_t_win_stop,f_start,f_stop,envelop_method,baseline_mode

#### basic imports
import sys,io,os,fnmatch,shutil

import matplotlib
matplotlib.use('PS')

#### nibabel import
#import nibabel as nib

##### nipype import
#from nipype import config
#config.enable_debug_mode()

import nipype
print nipype.__version__

import nipype.interfaces.io as nio

from nipype.interfaces.utility import IdentityInterface,Function
import nipype.pipeline.engine as pe

from neuropype_ephy.spectral import spectral_proc,plot_circular_connectivity

def get_freq_band(freq_band_name):

    from define_variables import freq_band_names,freq_bands
    
    if freq_band_name in freq_band_names:
        print freq_band_name
        print freq_band_names.index(freq_band_name)
        
        return freq_bands[freq_band_names.index(freq_band_name)]

####################################################################### by trigs ############################################################################################################

def create_infosource_merge_correct_cond_ampl_by_band():

    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id','trial_cond','freq_band_name','numPhase','trig_cond','exp']),name="infosource")
    
    #infosource.iterables = [('subject_id', subject_ids),('trial_cond',['all']),('freq_band_name',freq_band_names),('trig_cond',['odor','rest']),('exp',["E","R"])]

    ### test
    infosource.iterables = [('subject_id', ['VACJ']),('trial_cond',['all']),('freq_band_name',["beta"]),('numPhase',["1"]),('trig_cond',['odor']),('exp',["E"])]
        
    return infosource
    
def create_datasource_sigs_by_trigs():
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id','trial_cond','numPhase','trig_cond','exp'],outfields=['ts_file','channel_names_file','channel_coords_file']),name = 'datasource')
    datasource.inputs.base_directory = data_path
    datasource.inputs.template = '%s%s/%s/%s%s/%s%s%s%s%s'
    datasource.inputs.template_args = dict(
        ts_file=[["TS_" + exp ,"_all_cond_by_block_trigs",'subject_id','exp','numPhase',"correct_ts_",'trial_cond',"_by_",'trig_cond',"_trigs.npy"]],
        channel_names_file = [["TS_" + exp ,"_all_cond_by_block_trigs",'subject_id',"","","channel_names.txt","","","",""]],
        channel_coords_file = [["TS_" + exp ,"_all_cond_by_block_trigs",'subject_id',"","","channel_coords.txt","","","",""]]
        )
    datasource.inputs.sort_filelist = True
    return datasource
    
def create_main_workflow_spectral_modularity():

    main_workflow = pe.Workflow(name=correl_analysis_name)
    main_workflow.base_dir = main_path
    
    ## Info source
    infosource = create_infosource_merge_correct_cond_ampl_by_band()
        
    ## Data source
    datasource = create_datasource_sigs_by_trigs()

    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'trial_cond', datasource, 'trial_cond')
    main_workflow.connect(infosource, 'trig_cond', datasource, 'trig_cond')
    main_workflow.connect(infosource, 'numPhase', datasource, 'numPhase')
    main_workflow.connect(infosource, 'exp', datasource, 'exp')
    
    #### spectral
    spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],output_names = "conmat_file",function = spectral_proc),name = "spectral")
    spectral.inputs.con_method = con_method
    spectral.inputs.sfreq = sfreq
    
    #spectral.inputs.epoch_window_length = epoch_window_length
    main_workflow.connect(datasource, 'ts_file', spectral, 'ts_file')
    main_workflow.connect(infosource, ('freq_band_name',get_freq_band), spectral, 'freq_band')
    
    
    ##### plot spectral
    #plot_spectral = pe.MapNode(interface = Function(input_names = ["conmat_file","labels_file","nb_lines"],output_names = "plot_conmat_file",function = plot_circular_connectivity),iterfield = ['conmat_file'], name = "plot_spectral")
    
    ##plot_spectral.inputs.labels_file = MEG_elec_names_file
    #plot_spectral.inputs.nb_lines = 200
    
    #main_workflow.connect(spectral, "conmat_file", plot_spectral, 'conmat_file')
    #main_workflow.connect(datasource, "channel_names_file", plot_spectral, 'labels_file')
    
    
    
    
    ##### graph pipeline definition (density based) 
    #graph_den_pipe = create_pipeline_conmat_to_graph_density("graph_den_pipe",main_path,radatools_path,con_den = con_den, multi = False, mod = mod)
    
    #main_workflow.connect(spectral, "conmat_file", graph_den_pipe, 'compute_net_List.Z_cor_mat_file')
    
    #### only if "mod" is on 
    #if mod == True:
        #graph_den_pipe.inputs.community_rada.optim_seq = radatools_optim
    
        #graph_den_pipe.inputs.plot_igraph_modules_rada.coords_file = MEG_elec_coords_file 
        #graph_den_pipe.inputs.plot_igraph_modules_rada.labels_file = MEG_elec_names_file
            
    return main_workflow
    
if __name__ == '__main__':
    
    # run pipeline:
    main_workflow = create_main_workflow_spectral_modularity()
         
    ####### run
    
    #main_workflow.write_graph()
    main_workflow.config['execution'] = {'remove_unnecessary_outputs':'false'}
    main_workflow.run(plugin='MultiProc', plugin_args={'n_procs' : 8})    
    
    
    
