# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
import numpy as np

filename_results = 'All_cycles_features_{}.xls'

for exp in ['E']:
    writer = ExcelWriter(filename_results.format(exp))
    for subject in session.query(Subject).order_by('Subject.`group`'):
        for i,run in enumerate(subject.runs.filter_by(exp =exp).order_by('`index`')):
            print subject.name, exp

            if run.trials.count()==0:continue
            resp = run.respirationsignals[0]
            resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
            
            if resp.cycle_times is None:
                print 'compute resp'
                resp.cycle_times = compute_cycle_for_run(resp)
                resp.save(session)
    
            all_cycles = resp.cycle_times.magnitude
            n_cycles = all_cycles.shape[0]
            insp = all_cycles[:-1,0]
            expi = all_cycles[:-1,1]
            insp_next = all_cycles[1:,0]
            sr = resp.sampling_rate.rescale('Hz').magnitude
            
            insp_t, expi_t, insp_d, expi_d, insp_v, expi_v, insp_a, expi_a = [],[],[],[],[],[],[],[]
            for c in range(n_cycles-1):
                ix1 =int((insp[c]-resp.t_start.rescale('s').magnitude)*sr)
                ix2 =int((expi[c]-resp.t_start.rescale('s').magnitude)*sr)
                ix3 =int((insp_next[c]-resp.t_start.rescale('s').magnitude)*sr)
                insp_t.append(insp[c])
                expi_t.append(expi[c])
                insp_d.append(insp_next[c]-expi[c])
                expi_d.append(np.sum(resp_f[ix2:ix3])/sr)
                insp_v.append(np.max(resp_f[ix2:ix3]))
                expi_v.append(expi[c]-insp[c])
                insp_a.append(np.sum(resp_f[ix1:ix2])/sr)
                expi_a.append(np.min(resp_f[ix1:ix2]))

            df = DataFrame({})
            df['insp_time'] = insp_t
            df['exp_time'] = expi_t
            df['exp_duration'] = expi_d
            df['exp_volume'] = expi_v
            df['exp_max'] = expi_a
            df['insp_duration'] = insp_d
            df['insp_volume'] = insp_v
            df['insp_max'] = insp_a
            df.to_excel(writer, sheet_name=subject.name+'_E'+str(i+1))
    writer.save()