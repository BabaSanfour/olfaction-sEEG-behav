# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
from score_common import create_df_subject
from score_encoding import create_df_subject_encode
#from score_episodic import create_df_episodic_subject
#from score_reconnaissance import create_df_recognition_subject


import numpy as np

filename_results = 'Repiration features individual {}.xls'

time_windows = ['recognition', 'episodic', 'total' ]


filename_results = 'Repiration features individual {} {}.xls'

time_windowsE = ['total_encodage', ]
time_windowsR = ['recognition', 'episodic', 'total' ]


def get_respiration_feature_for_trial(trial, resp_f, t1, t2, n_cycle = 10):
    resp = trial.run.respirationsignals[0]
    if resp.cycle_times is None:
        print 'compute resp'
        resp.cycle_times = compute_cycle_for_run(resp)
        resp.save(session)
    
    cycles = resp.cycle_times.magnitude
    insp = cycles[:-1,0]
    expi = cycles[:-1,1]
    insp_next = cycles[1:,0]
    sr = resp.sampling_rate.rescale('Hz').magnitude
    
    f = features = OrderedDict()
    local_cycles,  = np.where((insp >=t1) & (insp<t2))
    
    
    list_mean = ['mean_insp_duration', 'mean_insp_volume', 'mean_insp_max',
                            'mean_exp_duration', 'mean_exp_volume', 'mean_exp_max',
                            'mean_cycle_duration',
                            ]
    for e in list_mean:
        f[e] = [ ]
    f['nb_cycle'] = local_cycles.size
    try:
        f['mean_frequency'] = local_cycles.size/(t2-t1)
    except:
        f['mean_frequency'] = None
    f['inv_mean_duration'] = None
    
    
    #~ for c, localcycle in enumerate(local_cycles[:n_cycle]):
    for c in range(n_cycle):
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
            
            
            f['mean_insp_duration'].append(f['{}_insp_duration'.format(c)])
            f['mean_insp_volume'].append(f['{}_insp_volume'.format(c)])
            f['mean_insp_max'].append(f['{}_insp_max'.format(c)])
            f['mean_exp_duration'].append(f['{}_exp_duration'.format(c)])
            f['mean_exp_volume'].append(f['{}_exp_volume'.format(c)])
            f['mean_exp_max'].append(f['{}_exp_max'.format(c)])
            f['mean_cycle_duration'].append(f['{}_insp_duration'.format(c)] + f['{}_exp_duration'.format(c)] )
        else :
            f['{}_insp_time'.format(c)] = np.nan
            f['{}_exp_time'.format(c)] = np.nan
            f['{}_exp_duration'.format(c)] = np.nan
            f['{}_exp_volume'.format(c)] = np.nan
            f['{}_exp_max'.format(c)] = np.nan
            f['{}_insp_duration'.format(c)] = np.nan
            f['{}_insp_volume'.format(c)] = np.nan
            f['{}_insp_max'.format(c)] = np.nan
    
    for e in list_mean:
        f[e] = np.mean(f[e])
    try:
        f['inv_mean_duration'] = 1./f['mean_cycle_duration']
    except:
        pass
    
    return features
    
    

###
def create_df_respiration_features_subject(subject, exp, time_window):
    prefix = time_window+'_'
    if exp =='R' and time_window == 'recognition':
        df = create_df_recognition_subject(subject)
    if exp == 'R' and time_window == 'episodic':
        df = create_df_episodic_subject(subject)
    if exp =='E':
        df = create_df_subject_encode(subject)
    else:
        df = create_df_subject(subject)
    ix = df.index
    for run in subject.runs.filter_by(exp =exp).order_by('`index`'):
        if run.trials.count()==0:continue
        resp = run.respirationsignals[0]
        resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
        for trial in run.trials:
            if time_window == 'recognition':
                t1 = trial.triggered_odor_time
                t2 = trial.recognition_time+trial.time
            elif time_window == 'episodic':
                t1 = trial.recognition_time+trial.time 
                t2 = trial.time+25.
            elif time_window == 'total':
                t1 = trial.triggered_odor_time
                t2 = trial.time+25.
            elif time_window == 'total_encodage':
                t1 = trial.triggered_odor_time
                #t2 = t1 + 4.5
                next_trial = run.trials.filter_by(index = trial.index+1).first()
                if next_trial is not None:
                    t2 = next_trial.time
                else:
                    t2 = trial.time
            features = get_respiration_feature_for_trial(trial, resp_f, t1, t2, n_cycle = 15)
            for k, v in features.items():
                if prefix+k not in df.columns:
                    df[prefix+k] = Series(index = ix)
                df.loc[trial.id,prefix+k] = v
    ## part added 01/06/2020
    cols = ['delta_time_next_trial','total_encodage_mean_exp_duration',
           'total_encodage_mean_exp_volume', 'total_encodage_mean_exp_max',
           'total_encodage_mean_insp_duration','total_encodage_mean_insp_volume',
           'total_encodage_mean_insp_max','trial_index']
    if exp == 'E':
        dict_agg = {}
        for key in cols:
            dict_agg[key] = 'mean' if key != 'trial_index' else 'count'
        df = df.groupby(['subject_name','odor_num']).agg(dict_agg)
    elif exp =='R':
        cols2 = [ c.replace('encodage_','') for c in cols[1:-1]]
        dict_agg = {}
        for key in cols2:
            dict_agg[key] = 'mean' 
        df = df.groupby(['subject_name','odor_num']).agg(dict_agg)
    #~ df = df[df.odor_num!=17]
    return df



def test1():
    subject = session.query(Subject)[5]
    for time_window in time_windows:
        df = create_df_respiration_features_subject(subject, time_window)
        print df

def export_xl_respiration_features():

    exps, time_windows = ['R'], ['total']
    for exp, time_window in zip(exps,time_windows):
        print(exp,time_window)
        writer = ExcelWriter(filename_results.format(time_window, exp))
        for subject in session.query(Subject).order_by('Subject.`group`'):
            print subject.name, exp
            df = create_df_respiration_features_subject(subject,exp, time_window)
            df.to_excel(writer, sheet_name=subject.name)
        writer.save()    






if __name__ == '__main__':
    #~ test1()
    
    export_xl_respiration_features()
    



