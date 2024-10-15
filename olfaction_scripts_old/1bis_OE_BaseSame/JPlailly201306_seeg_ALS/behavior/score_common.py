# -*- coding: utf-8 -*-
from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

import codes

from collections import OrderedDict
from respiration_tools import compute_cycle_for_run, get_time_first_inspi_or_trig

    
def create_df_subject(subject):
    
    ix = [ ]
    for run in subject.runs.filter_by(exp = 'R'):
        for trial in run.trials:
            ix.append(trial.id)
    n = len(ix)
    df = DataFrame(index = ix)
    headers =  ['subject_name','subject_num', 'subject_group','run_index', 'trial_index',
                                'trial_time', 'trial_triggered_odor_time', 'trial_recognition', 
                                'trial_context_time', 'trial_click_image', 'trial_click_posistion',
                                'odor_num', 'odor_name', 'odor_pleasantness', 
                                'image_context', 'context_day',
                                'is_target' , 'repetition_of_target', 'repetition_of_distractor',
                                'first_inspi_or_trig',
                                'trigger_is_during_inspi', 'trigger_is_first_part_inspi',
                                'delta_time_next_trial', 'delta_time_clic_expi',
                                'first_expi_duration','duration_odor_stim','time_before_trig_odor',
                                ]

    for h in headers:
        if h in ['subject_name', 'subject_num', 'subject_group', 'odor_name',  'image_context']:
            df[h] = Series(index = ix, dtype = 'S')
        elif h in ['repetition_of_target', 'repetition_of_distractor', 'trial_click_posistion']:
            df[h] = Series(index = ix, dtype = 'i')
        else:
            df[h] = Series(index = ix)
    
    df['subject_name'] = Series([str(subject.name)]*n, index = ix)
    df['subject_num'] = Series([str(subject.num)]*n, index = ix)
    df['subject_group'] = Series([str(subject.group)]*n, index = ix)

    for run in subject.runs.filter_by(exp = 'R'):
        resp = run.respirationsignals[0]
        if resp.cycle_times is None:
            print 'compute resp'
            resp.cycle_times = compute_cycle_for_run(resp)
            resp.save(session)
        cycles = resp.cycle_times.magnitude
        insp = cycles[:-1,0]
        expi = cycles[:-1,1]
        
        for trial in run.trials:
            df.loc[trial.id, 'run_index'] = run.index
            df.loc[trial.id, 'trial_index'] = trial.index
            df.loc[trial.id, 'trial_time'] = trial.time
            df.loc[trial.id, 'trial_triggered_odor_time'] = trial.triggered_odor_time
            df.loc[trial.id, 'trial_context_time'] = trial.context_time
            df.loc[trial.id, 'trial_recognition'] = trial.recognition
            df.loc[trial.id, 'trial_click_image']= trial.selected_image
            df.loc[trial.id, 'trial_click_posistion'] = trial.click_posistion
            
            
            df.loc[trial.id, 'odor_num'] = '{:03d}'.format( trial.odor_num )
            df.loc[trial.id, 'odor_name'] = unicode(trial.odor_name)
            if subject.odor_pleasantness is not None:
                df.loc[trial.id, 'odor_pleasantness'] = subject.odor_pleasantness[trial.odor_num-1] - 5.
            
            if trial.odor_num in codes.num_odor_to_ref[subject.name]:
                df.loc[trial.id, 'image_context'] = str(codes.num_odor_to_ref[subject.name][trial.odor_num][3])
                df.loc[trial.id, 'is_target'] =  1
                df.loc[trial.id, 'context_day'] = subject.encoding_order.index(df['image_context'][trial.id])+1
            else:
                df.loc[trial.id, 'image_context'] = None
                df.loc[trial.id, 'is_target'] =  0
            
            #~ if (trial.run.group == '1' and trial.odor_num<10) or (trial.run.group == '2' and trial.odor_num>=10):
                #~ df['is_target'][trial.id] =  1
                #~ df['context_day'][trial.id] = subject.encoding_order.index(df['image_context'][trial.id])+1
            #~ elif (trial.run.group == '1' and trial.odor_num>=10) or trial.run.group == '2' and trial.odor_num<10:
                #~ df['is_target'][trial.id] =  0
            
            # is trigger OK or not
            df.loc[trial.id, 'first_inspi_or_trig'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle)
            trial.first_inspiration_time = df.loc[trial.id, 'first_inspi_or_trig']
            t = trial.triggered_odor_time
            ind, = where((t>insp) & (t<=expi))
            if ind.size == 1:
                percent = (t-insp[ind[0]])/(expi[ind[0]]-insp[ind[0]])
                if percent <max_percent_for_respiration_cycle:
                    df.loc[trial.id, 'trigger_is_during_inspi'] = True
                    df.loc[trial.id, 'trigger_is_first_part_inspi'] = True
                else:
                    df.loc[trial.id, 'trigger_is_during_inspi'] = True
                    df.loc[trial.id, 'trigger_is_first_part_inspi'] = False
            else:
                    df.loc[trial.id, 'trigger_is_during_inspi'] = False
                    df.loc[trial.id, 'trigger_is_first_part_inspi'] = False

            df.loc[trial.id, 'delta_time_clic_expi'] = trial.triggered_odor_time-trial.time
            #print trial.triggered_odor_time, trial.time, trial.first_inspiration_time
            df.loc[trial.id,'first_expi_duration'] = trial.first_inspiration_time - trial.triggered_odor_time
            df.loc[trial.id,'duration_odor_stim'] = 4.5- (trial.first_inspiration_time - trial.triggered_odor_time)
            df.loc[trial.id,'time_before_trig_odor'] = trial.first_inspiration_time - trial.time
            
                
    for num in np.unique(df[df.is_target==1]['odor_num']):
        sel = df['odor_num'] == num
        df.loc[sel, 'repetition_of_target'] = np.arange(sum(sel))
    for num in np.unique(df[df.is_target==0]['odor_num']):
        sel = df['odor_num'] == num
        df.loc[sel, 'repetition_of_distractor'] = np.arange(sum(sel))
    
            

    session.commit()
    
    return df
    
def test_create_df_subject():
    #~ subject = session.query(Subject)[1]
    subject = session.query(Subject).filter_by(name = 'LEFC').first()
    df =  create_df_subject(subject)
    #~ print df[ ['is_target', 'image_context']]
    print df



    
if __name__ =='__main__':
    test_create_df_subject()


