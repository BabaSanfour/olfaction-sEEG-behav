# -*- coding: utf-8 -*-
from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

import codes
from collections import OrderedDict
from respiration_tools import compute_cycle_for_run, get_time_first_inspi_or_trig

pos_files_dirname = './FichiersPOS'
if not os.path.exists(pos_files_dirname):
    os.mkdir(pos_files_dirname)

filename_results1 = 'Fichier_POS_Encoding.xls'

def create_df_encoding_pos(subject):
    
    ix = [ ]
    for run in subject.runs.filter_by(exp='E').order_by('run_index'):
        for trial in run.trials:
            ix.append(trial.id)
    n = len(ix)
    df = DataFrame(index = ix)
    
    df['encoding_day'] = Series(index = ix, dtype = 'S')
    for k in ['start', 'odor', 'first_inspi']:
        df[k] = Series(index = ix, dtype = int)
        df[k+'_code'] = Series(index = ix, dtype = int)

    for run in subject.runs.filter_by(exp='E'):
        if run.trials.count()==0:continue
        resp = run.respirationsignals[0]
        if resp.cycle_times is None:
            print 'compute resp'
            resp.cycle_times = compute_cycle_for_run(resp)
            resp.save(session)
        cycles = resp.cycle_times.magnitude
        insp = cycles[:-1,0]
        expi = cycles[:-1,1]
        
        for trial in run.trials:
            df.loc[trial.id,'encoding_day'] = run.index
  
            df.loc[trial.id,'start'] = trial.time*512.
            df.loc[trial.id,'start_code'] = trial.odor_num+100
            
            df.loc[trial.id,'odor'] = trial.triggered_odor_time*512.
            df.loc[trial.id,'odor_code'] = trial.odor_num+200
            
            df.loc[trial.id,'first_inspi'] = trial.first_inspiration_time*512.
            df.loc[trial.id,'first_inspi_code'] = trial.odor_num+300
            
    session.commit()
    print df
    return df

def export_encoding_pos_all():
    
    writer = ExcelWriter('Fichier_POS_Encoding.xls')
    for subject in session.query(Subject).order_by('Subject.`group`'):
        print subject.name
        df = create_df_encoding_pos(subject)
        df.to_excel(writer, sheet_name=subject.name)
        
        for encoding_day in [1,2]:
            all = []
            for k in ['start', 'odor', 'first_inspi', ]:
                df2 = df[df['encoding_day']==encoding_day]
                df2 = df2[[k, k+'_code']]
                df2.columns = ['sample', 'code']
                df2['reject'] = 0
                all.append(df2)
            all = concat(all, axis=0)
            all = all[~all['sample'].isnull()]
            all.to_csv(os.path.join(pos_files_dirname,'{}_{}{}.pos'.format(subject.name, 'E', encoding_day)), header = False, index = False, sep = '\t', float_format = '%.0f' )
        
        
        
    writer.save()    
    
    
if __name__ =='__main__':
    export_encoding_pos_all()
    



