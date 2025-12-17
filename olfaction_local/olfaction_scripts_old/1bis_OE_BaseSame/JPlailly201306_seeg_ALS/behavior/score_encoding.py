# -*- coding: utf-8 -*-



from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
#from score_reconnaissance import get_all as get_all_rec
#from score_episodic import get_all as get_all_epi


import codes

from collections import OrderedDict
from respiration_tools import compute_cycle_for_run, get_time_first_inspi_or_trig



filename_results1 = 'encoding_individual_results.xls'
filename_results2 = 'encoding_summary.xls'


show_figures = False

def create_df_subject_encode(subject):
    
    ix = [ ]
    for run in subject.runs.filter_by(exp='E').order_by('run_index'):
        for trial in run.trials:
            ix.append(trial.id)
    n = len(ix)
    df = DataFrame(index = ix)
    headers =  ['subject_name', 'subject_group','run_index',  'trial_index',
                                'trial_time', 
                                'odor_num', 'odor_name', 'odor_pleasantness',
                                'delta_time_next_trial', 'delta_time_clic_expi',
                                'image_context', 'encoding_day','trial_triggered_odor_time',
                                'first_inspi_or_trig','time_no_odor_btw_trials','first_expi_duration',
                                'duration_odor_stim', 'time_before_trig_odor',
                                'trigger_is_during_inspi', 'trigger_is_first_part_inspi',
                                
                                ]

    for h in headers:
        if h in ['subject_name', 'subject_num', 'subject_group', 'odor_name',  'image_context']:
            df[h] = Series(index = ix, dtype = 'S')
        elif h in [ ]:
            df[h] = Series(index = ix, dtype = 'i')
        else:
            df[h] = Series(index = ix)
    
    df['subject_name'] = Series([str(subject.name)]*n, index = ix)
    df['subject_group'] = Series([str(subject.group)]*n, index = ix)

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
            df.loc[trial.id, 'run_index'] = run.index
            df.loc[trial.id,'trial_index'] = trial.index
            df.loc[trial.id,'trial_time'] = trial.time
            
            
            df.loc[trial.id,'odor_num'] = '{:03d}'.format( trial.odor_num )
            df.loc[trial.id,'odor_name'] = unicode(trial.odor_name)
            
            
            if subject.odor_pleasantness is not None:
                df['odor_pleasantness'][trial.id] = subject.odor_pleasantness[trial.odor_num-1]
            
            next_trial = run.trials.filter_by(index = trial.index+1).first()
            if next_trial is not None:
                df.loc[trial.id,'delta_time_next_trial'] = next_trial.time - trial.time
                df.loc[trial.id,'time_no_odor_btw_trials'] = next_trial.time - (trial.triggered_odor_time+4.5)
            
            df.loc[trial.id,'image_context'] = run.image
            
            
            df.loc[trial.id,'encoding_day'] = run.index
            df.loc[trial.id,'trial_triggered_odor_time'] = trial.triggered_odor_time
            df.loc[trial.id, 'delta_time_clic_expi'] = trial.triggered_odor_time-trial.time
            
            # is trigger OK or not
            df.loc[trial.id,'first_inspi_or_trig'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle)
            df.loc[trial.id,'first_expi_duration'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle) - trial.triggered_odor_time
            df.loc[trial.id,'duration_odor_stim'] = 4.5-(get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle) - trial.triggered_odor_time)
            df.loc[trial.id,'time_before_trig_odor'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle)-(trial.time)

            trial.first_inspiration_time = df.loc[trial.id,'first_inspi_or_trig']
    
            t = trial.triggered_odor_time
            ind, = where((t>insp) & (t<=expi))
            if ind.size == 1:
                percent = (t-insp[ind[0]])/(expi[ind[0]]-insp[ind[0]])
                if percent <max_percent_for_respiration_cycle:
                    df.loc[trial.id,'trigger_is_during_inspi'] = True
                    df.loc[trial.id,'trigger_is_first_part_inspi'] = True
                else:
                    df.loc[trial.id,'trigger_is_during_inspi'] = True
                    df.loc[trial.id,'trigger_is_first_part_inspi'] = False
            else:
                    df.loc[trial.id,'trigger_is_during_inspi'] = False
                    df.loc[trial.id,'trigger_is_first_part_inspi'] = False
                

    session.commit()

    return df[df > 0]

    
def test_create_df_subject_encode():
    subject = session.query(Subject)[5]
    print subject
    df =  create_df_subject_encode(subject)
    print df
    #~ print df[ ]


def export_xl_encoding_scores_individual():
    from respiration_features import create_df_respiration_features_subject
    writer = ExcelWriter(filename_results1)
    for subject in session.query(Subject).order_by('Subject.`group`'):
        print subject.name
        df = create_df_respiration_features_subject(subject,'E', 'total_encodage')
        df.to_excel(writer, sheet_name=subject.name)
        #~ df2 = df[['subject_name', 'odor_num', 'score_recognition']]
        #~ df2.to_excel(writer, sheet_name=subject.name)
    writer.save()    
    
def get_all():
    from respiration_features import create_df_respiration_features_subject
    all = [ ] 
    for subject in session.query(Subject).order_by('Subject.`group`'):
        df = create_df_respiration_features_subject(subject, 'E', 'total_encodage')
        #print df.columns
        all.append(df)
    return concat(all)
    

def export_xl_encoding_scores_synthesis():
    writer = ExcelWriter(filename_results2)
    
    #Encoding
    all = get_all()
    # par sujet
    select = all.groupby(['subject_name','odor_num']).mean()#
    #~ select = all.groupby(['subject_name', 'image_context', 'encoding_day'])
    
    # l = [
    #         #~ select['image_context'].first(),
    #         #~ select['encoding_day'].first(),
    #         select['trial_index'].count(),
    #         select['delta_time_next_trial'].mean(),
    #         select['first_expi_duration'].mean(),
    #         select['first_expi_duration'].min(),
    #         select['time_no_odor_btw_trials'].min(),
    #         select['duration_odor_stim'].mean(),
    #         select['duration_odor_stim'].min(),
    #         select['duration_odor_stim'].max(),
    #         select['duration_odor_stim'].std(),
    #         select['time_before_trig_odor'].mean(),
    #         select['time_before_trig_odor'].min(),
    #         select['time_before_trig_odor'].max(),
    #         select['time_before_trig_odor'].std(),
    #         select['delta_time_clic_expi'].mean(),
    #         select['delta_time_clic_expi'].min(),
    #         select['delta_time_clic_expi'].max(),
    #         select['delta_time_clic_expi'].std()
    #         ]
    
    # resp_features = ['mean_insp_duration', 'mean_insp_volume', 'mean_insp_max',
    #                         'mean_exp_duration', 'mean_exp_volume', 'mean_exp_max',
    #                         'mean_cycle_duration','nb_cycle', 'mean_frequency',
    #                         'inv_mean_duration', 
    #                         ]
    # for e in resp_features:
    #     l.append(select['total_encodage_'+e].mean())
        
    #df_encode = concat(l, axis = 1)
    #~ df_encode.columns = ['image_context','encoding_day', 'count', 'delta_time_next_trial' ] + resp_features
    # df_encode.columns = ['count', 'delta_time_next_trial','first_expi_duration','first_expi_duration_min',
    #                     'time_no_odor_min','odor_stim_mean','odor_stim_min','odor_stim_max','odor_stim_std',
    #                     'before_odor_mean','before_odor_min','before_odor_max','before_odor_std',
    #                     'clic_expi_mean','clic_expi_min','clic_expi_max','clic_expi_std'] + resp_features
    # #~ df_encode= select.agg([ 'count'])#['image_context', 'encoding_day']
    #~ df_encode= select.value_counts()
    #~ df_encode = DataFrame(select.size())
    #~ print df_encode
    
    #~ #Rappel Recognition
    #~ all2 = get_all_rec()
    #~ select = all2.groupby(['subject_name', 'odor_num' ])
    #~ df_rec = DataFrame(select['score_recognition'].first())
    #~ df_rec.columns = ['score_recognition']
    
    #~ #Rappel Episodic
    #~ all3 = get_all_epi()
    #~ select = all3.groupby(['subject_name', 'odor_num' ])
    #~ df_epi = DataFrame(select['score_episodic_strict'].first())
    #~ df_epi.columns = ['score_episodic_strict']
    
    #~ df_encode.to_excel(writer, sheet_name='df_encode')
    #~ df_rec.to_excel(writer, sheet_name='df_rec')
    #~ df_epi.to_excel(writer, sheet_name='df_epi')
    #~ df = concat([ df_encode, df_rec, df_epi], axis = 1, join = 'inner')
    #~ print df
    
    #~ df.to_excel(writer, sheet_name = 'encoding_rappel')
    select.to_excel(writer, sheet_name = 'encoding_rappel')
    writer.save()
    
    if show_figures:
        pass
        #~ plot_df(score_by_subject, score_episodic_strict_list)
        #~ plot_df(score_by_subject2, score_episodic_elargi_list)
        #~ plot_df(score_by_odor1, score_episodic_strict_list)
        #~ plot_df(score_by_odor2, score_episodic_elargi_list)
        #~ plot_df(score_by_context_day1, score_episodic_strict_list)
        #~ plot_df(score_by_context_day2, score_episodic_elargi_list)
        #~ plot_df(score_by_image_context1, score_episodic_strict_list)
        #~ plot_df(score_by_image_context2, score_episodic_elargi_list)
        
        #~ pyplot.show()
    
if __name__ =='__main__':
    #~ test_create_df_subject_encode()
    #export_xl_encoding_scores_individual()
    
    export_xl_encoding_scores_synthesis()



