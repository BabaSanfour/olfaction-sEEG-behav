# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
from score_common import create_df_subject
from score_encoding import create_df_subject_encode
from itertools import product

import numpy as np

path2save = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/1bis_OE_BaseSam/JPlailly201306_seeg_ALS/behavior/respiration_amplitude/'
filename_results = 'Respi_ind_{}_{}.xls'
time_windows = ['total_encodage']
exps = ['E']

def get_respiration_feature_for_trial(trial, resp_f, t1, t2, n_cycle = 2):
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
    local_cycles, = np.where((insp >=t1) & (insp<t2))
    f['nb_cycle'] = local_cycles.size
  
    for c in range(n_cycle):
        if c < local_cycles.size:
            localcycle = local_cycles[c]
            ix1 =int((insp[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            ix2 =int((expi[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            ix3 =int((insp_next[localcycle]-resp.t_start.rescale('s').magnitude)*sr)
            f['{}_insp_time'.format(c)] = insp[localcycle]
            f['{}_exp_time'.format(c)] = expi[localcycle]
            f['{}_exp_ampl'.format(c)] = np.max(resp_f[ix2:ix3])
            f['{}_insp_max'.format(c)] = np.max(resp_f[ix1:ix2])
            
        else :
            f['{}_insp_time'.format(c)] = np.nan
            f['{}_exp_time'.format(c)] = np.nan
            f['{}_exp_ampl'.format(c)] = np.nan
            f['{}_insp_max'.format(c)] = np.nan
    

    # try:
    #     f['inv_mean_duration'] = 1./f['mean_cycle_duration']
    # except:
    #     pass
    
    return features

###
def create_df_respiration_features_subject(subject, exp, time_window):
    prefix = time_window+'_'
    if exp =='E':
        df = create_df_subject_encode(subject)
    ix = df.index
    for run in subject.runs.filter_by(exp =exp).order_by('`index`'):
        if run.trials.count()==0:continue
        resp = run.respirationsignals[0]
        resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
        for trial in run.trials:
            if time_window == 'total_encodage':
                t1 = trial.triggered_odor_time
                t2 = t1 + 4
            features = get_respiration_feature_for_trial(trial, resp_f, t1, t2, n_cycle = 2)
            for k, v in features.items():
                if prefix+k not in df.columns:
                    df[prefix+k] = Series(index = ix)
                df.loc[trial.id,prefix+k] = v
    #~ df = df[df.odor_num!=17]
    return df

def export_xl_respiration_features():
    for exp, time_window in product(exps, time_windows):
        writer = ExcelWriter(path2save+filename_results.format(time_window, exp))
        for subject in session.query(Subject).order_by('Subject.`group`'):
            print subject.name, exp
            df = create_df_respiration_features_subject(subject,exp,  time_window)
            df.to_excel(writer, sheet_name=subject.name)
        writer.save()    


if __name__ == '__main__':
    export_xl_respiration_features()