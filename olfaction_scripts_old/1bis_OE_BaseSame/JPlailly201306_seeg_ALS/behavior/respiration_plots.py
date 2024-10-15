# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
from score_common import create_df_subject
from score_encoding import create_df_subject_encode
import numpy as np

for subject in session.query(Subject).order_by('Subject.`group`'):
    times_all, resp_all = np.array([]),np.array([])
    print subject.name
    for run in subject.runs.filter_by(exp ='E').order_by('`index`'):
        resp = run.respirationsignals[0]
        sig = resp.signal
        t_start = resp.t_start.rescale('s').magnitude
        sr = resp.sampling_rate.rescale('Hz').magnitude
        times = np.arange(sig.shape[0])/sr + t_start
        for trial in run.trials:
            t1 = trial.triggered_odor_time
            t2 = t1 + 4.5
            sig_chunk = sig[(times>=t1) & (times<=t2)][np.newaxis]
            resp_all = np.vstack((resp_all,sig_chunk)) if np.size(resp_all) else sig_chunk
            print(sig_chunk.shape,resp_all.shape)
            times_chunk = times[(times>=t1) & (times<=t2)]-t1
            times_chunk = times_chunk[np.newaxis]
            times_all = np.vstack((times_all,times_chunk)) if np.size(times_all) else times_chunk
            print(times_chunk.shape,times_all.shape)
        print('run',resp_all.shape, times_all.shape)
    print('all trials',resp_all.shape, times_all.shape)