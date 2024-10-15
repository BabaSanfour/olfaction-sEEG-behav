# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
from score_common import create_df_subject
from score_encoding import create_df_subject_encode
from score_episodic import create_df_episodic_subject
from score_reconnaissance import create_df_recognition_subject

import numpy as np

filename_results = 'Repiration features individual {} {}.xls'
time_windowsE = ['odor1_encoding','odor2_encoding','no_odor_encoding']

def get_respiration_feature_for_trial(trial, resp_f, t0, t1, t2):
    print t0,t1,t2
    resp = trial.run.respirationsignals[0]
    if resp.cycle_times is None:
        print 'compute resp'
        resp.cycle_times = compute_cycle_for_run(resp)
        resp.save(session)
    
    all_cycles = resp.cycle_times.magnitude
    n_cycles = all_cycles.shape[0]
    insp = all_cycles[:-1,0]
    expi = all_cycles[:-1,1]
    insp_before = np.insert(all_cycles[0:-2,0],0,0)
    insp_next = all_cycles[1:,0]
    insp_next2 = all_cycles[2:,0]
    insp_next2 = np.insert(insp_next2,-1,insp_next2[-1]+10)
    print insp[:10],insp_next2[:10]
    sr = resp.sampling_rate.rescale('Hz').magnitude
    
    f = features = OrderedDict()
    no_odor_cycles, = np.where((insp<t0)&(insp>insp_before)) #cycles between odor periods
    trigg_cycles,  = np.where((expi<=t1)&(t1<expi+0.2)) #first cycle triggering odor
    odor0_cycles,  = np.where(insp>t1)[0]#&(insp_next2<t1)) #first cycle while odor is sent
    print no_odor_cycles,trigg_cycles,odor0_cycles
    #for c, localcycle in enumerate(local_cycles[:n_cycle]):

    for c in range(n_cycles):
        if c < local_cycles.size:
            localcycle = local_cycles[c]
            ix1 =int((insp[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            ix2 =int((expi[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            ix3 =int((insp_next[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            f['{}_insp_time'.format(c)] = insp[localcycle]
            f['{}_exp_time'.format(c)] = expi[localcycle]
            f['{}_exp_duration'.format(c)] = insp_next[localcycle]-expi[localcycle]
            f['{}_exp_volume'.format(c)] = np.sum(resp_f[ix2:ix3])/sr
            f['{}_exp_max'.format(c)] = np.max(resp_f[ix2:ix3])
            f['{}_insp_duration'.format(c)] = expi[localcycle]-insp[localcycle]
            f['{}_insp_volume'.format(c)] = np.sum(resp_f[ix1:ix2])/sr
            f['{}_insp_max'.format(c)] = np.min(resp_f[ix1:ix2])
       
    return features
   

###
def create_df_respiration_features_subject(subject, exp, time_window):
    prefix = time_window+'_'
    print exp, time_window
    df = create_df_subject_encode(subject)
    ix = df.index
    for run in subject.runs.filter_by(exp =exp).order_by('`index`'):
        if run.trials.count()==0:continue
        resp = run.respirationsignals[0]
        resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
        for trial in run.trials:
            t0 = trial.time
            t1 = trial.triggered_odor_time
            t2 = trial.triggered_odor_time + 4.
            features = get_respiration_feature_for_trial(trial,resp_f,t0,t1,t2)
    #         for k, v in features.items():
    #             if prefix+k not in df.columns:
    #                 df[prefix+k] = Series(index = ix)
    #             df.loc[trial.id,prefix+k] = v
    # return df

def export_xl_respiration_features():

    for exp, time_window in zip(['E'],['odor1_encoding']):
        writer = ExcelWriter(filename_results.format(time_window, exp))
        for subject in session.query(Subject).order_by('Subject.`group`'):
            print subject.name, exp, time_window
            df = create_df_respiration_features_subject(subject, exp, time_window)
            df.to_excel(writer, sheet_name=subject.name)
        writer.save()    
    





if __name__ == '__main__':
    
    export_xl_respiration_features()
    



