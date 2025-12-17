# -*- coding: utf-8 -*-

from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

from respiration_tools import get_time_first_inspi_or_trig


from score_common import create_df_subject

from collections import OrderedDict

from matplotlib import pyplot



score_recognition_list = [ 'hit',  'miss', 'cr', 'fa' ]
score_colors = {'hit' :'#0080FF' ,
                                'miss':'#CEE3F6',
                                'cr':'#40FF00', 
                                'fa' :'#FF8000',  
                                }

filename_results1 = 'recognition_individual_results .xls'
filename_results2 = 'recognition_summary.xls'
filename_results3 = 'recognition_timing.xls'



#~ show_figures = True
show_figures = False

	

def score(ref, ans):
    if ans is None:
        return None
    elif (ref == 1) and (ans ==1):
        return 'hit'
    elif (ref == 0) and (ans ==1):
        return 'fa'
    elif (ref == 1) and (ans ==0):
        return 'miss'
    elif (ref == 0) and (ans ==0):
        return 'cr'
    
    
def create_df_recognition_subject(subject):
    df =  create_df_subject(subject)
    ix = df.index
    df['recognition_time'] = Series(index = ix)
    df['score_recognition']  = Series(index = ix, dtype = 'S')
    df['first_inspi_or_trig']  = Series(index = ix)
    df['recognition_rt']  = Series(index = ix)
    for s in score_recognition_list:
        df[s] =  Series(index = ix)
    for run in subject.runs.filter_by(exp = 'R'):
        resp = run.respirationsignals[0]
        for trial in run.trials:
            df.loc[trial.id,'recognition_time'] = trial.recognition_time
            trial.score_recognition = score(df.loc[trial.id, 'is_target'], trial.recognition)
            df.loc[trial.id,'score_recognition'] =  u'' if trial.score_recognition is None else trial.score_recognition
            df.loc[trial.id, 'first_inspi_or_trig'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle)
            trial.first_inspiration_time = df.loc[trial.id, 'first_inspi_or_trig']
            df.loc[trial.id, 'recognition_rt'] = trial.recognition_time - (trial.triggered_odor_time-trial.time)
            df.loc[trial.id,'recognition_rt'] = trial.recognition_time - (df.loc[trial.id, 'first_inspi_or_trig'] -trial.time)
            
            for s in score_recognition_list:
                df.loc[trial.id,s] = int(s == trial.score_recognition)

    session.commit()
    return df
    
def test_create_df_recognition_subject():
    #~ subject = session.query(Subject)[5]
    subject = session.query(Subject).filter_by(name = 'LEFC').first()
    df =  create_df_recognition_subject(subject)
    print df
    print df[ ['is_target', 'score_recognition', 'image_context']]
    print df[score_recognition_list]

def export_xl_recognition_scores_individual():
    from respiration_features import create_df_respiration_features_subject
    writer = ExcelWriter(filename_results1)
    for subject in session.query(Subject).order_by('Subject.`group`'):
        print subject.name
        df =  create_df_respiration_features_subject(subject,'R','total')
        df.to_excel(writer, sheet_name=subject.name)
    writer.save()    

def get_all():
    from respiration_features import create_df_respiration_features_subject
    all = [ ]  
    for subject in session.query(Subject).order_by('Subject.`group`'):
        print subject.name
        df =  create_df_respiration_features_subject(subject,'R','total') 
        all.append(df)
    return concat(all)

def test_get_all():
    all = get_all()
    print all
    
def add_tr_by_score(all, select, df):
    for s in score_recognition_list:
        df[s+'_rt'] = Series(index = df.index)
    for name, group in select:
        for s in score_recognition_list:
            sub = group[group['score_recognition'] == s]
            df[s+'_rt'][name] = sub['recognition_rt'].mean()

def export_xl_timing_recognition_synthesis():
    writer = ExcelWriter(filename_results3)

    all = get_all()
    print all.columns
    all = all[all['is_target'] ==1]
    select = all.groupby(['subject_name'])#
    #~ select = all.groupby(['subject_name', 'image_context', 'encoding_day'])
    
    l = [
            select['trial_index'].count(),
            select['delta_time_next_trial'].mean(),
            select['first_expi_duration'].mean(),
            select['duration_odor_stim'].mean(),
            select['duration_odor_stim'].min(),
            select['duration_odor_stim'].max(),
            select['duration_odor_stim'].std(),
            select['time_before_trig_odor'].mean(),
            select['time_before_trig_odor'].min(),
            select['time_before_trig_odor'].max(),
            select['time_before_trig_odor'].std(),
            select['delta_time_clic_expi'].mean(),
            select['delta_time_clic_expi'].min(),
            select['delta_time_clic_expi'].max(),
            select['delta_time_clic_expi'].std(),
            select['recognition_rt'].mean()
            ]
    
    resp_features = ['mean_insp_duration', 'mean_insp_volume', 'mean_insp_max',
                            'mean_exp_duration', 'mean_exp_volume', 'mean_exp_max',
                            'mean_cycle_duration','nb_cycle', 'mean_frequency',
                            'inv_mean_duration', 
                            ]
    for e in resp_features:
        l.append(select['recognition_'+e].mean())
        
    df_encode = concat(l, axis = 1)
    #~ df_encode.columns = ['image_context','encoding_day', 'count', 'delta_time_next_trial' ] + resp_features
    df_encode.columns = ['count', 'delta_time_next_trial','first_expi_duration',
                        'odor_stim_mean','odor_stim_min','odor_stim_max','odor_stim_std',
                        'before_odor_mean','before_odor_min','before_odor_max','before_odor_std',
                        'clic_expi_mean','clic_expi_min','clic_expi_max','clic_expi_std','tr_rec'] + resp_features
    df_encode.to_excel(writer, sheet_name = 'recognition_timing')
    writer.save()
    


def export_xl_recognition_scores_synthesis():
    all = get_all()
    print all.columns
    writer = ExcelWriter(filename_results2)
    
    # #par sujet
    # all2 = all[all['is_target'] ==1]
    # select = all2.groupby(['subject_name',])
    # score_by_subject = select.sum() [score_recognition_list]

    # add_tr_by_score (all, select, score_by_subject)
    # score_by_subject.to_excel(writer, sheet_name='by_subject')
    # print score_by_subject[score_recognition_list]

    # # par odeur
    # all2 = all[all['is_target'] ==1]
    # select = all2.groupby(['subject_name','odor_num'])
    # score_by_subject = select.sum() [score_recognition_list]

    # add_tr_by_score (all, select, score_by_subject)
    # score_by_subject.to_excel(writer, sheet_name='by_odor')
    # print score_by_subject[score_recognition_list]

    # par odeur - JUST FOR BREATHING MODULATIONS
    all2 = all[all['is_target'] ==1]
    select = all2.groupby(['subject_name','odor_num']).mean()
    select.to_excel(writer, sheet_name='by_odor')

    # # par context_day
    # all2 = all[all['is_target'] ==1]
    # select = all2.groupby(['subject_name', 'context_day'])
    # score_by_context_day = select.sum()[['hit', 'miss']]
    # add_tr_by_score (all, select, score_by_context_day)
    # score_by_context_day.to_excel(writer, sheet_name='by_context_day')
    # print  score_by_context_day[['hit', 'miss']]
    
    # # par image context
    # all2 = all[all['is_target'] ==1]
    # select = all2.groupby(['subject_name', 'image_context'])
    # score_by_image_context = select.sum()[['hit', 'miss']]
    # add_tr_by_score (all, select, score_by_image_context)
    # score_by_image_context.to_excel(writer, sheet_name='by_image_context')
    # print score_by_image_context[['hit', 'miss']]
    
    writer.save()
    
    if show_figures:
        #~ plot_score_df(score_by_subject, title = 'sujet')
        #~ plot_score_df(score_by_odor, title = 'sujet/odor')
        #~ plot_score_df(score_by_repetition_of_target, title = 'sujet/repetition')
        #~ plot_score_df(score_by_repetition_of_distractor)
        #~ plot_score_df(score_by_context_day, title = 'sujet/context_day')
        #~ plot_score_df(score_by_image_context, title = 'sujet/image_context')

        pyplot.show()

def plot_score_df(df, title = ''):
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title(title)
    
    l = []
    for i,index in enumerate(df.index):
        if type(index) == tuple:
            l.append(('{} '*len(index)).format(*index))
        else:
            l.append(index)
        bottom = 0
        for s in score_recognition_list:
            if s in df:
                y = df[s][index]
                ax.bar(i, y, width =.8, color = score_colors[s], bottom = bottom)
                if y>0:
                    ax.text(i+.4, bottom+y/2, s,horizontalalignment='center',verticalalignment='center')
                bottom += y
    ax.set_xticks(arange(len(l))+.4)
    ax.set_xticklabels(l, rotation = 45)


def plot_confidence_df(df, suffix = '_confidence_avg', title = '', same_axe = False):
    fig = pyplot.figure()
    
    if same_axe:
        ax = fig.add_subplot(1,1,1)
    
    sp = 1
    for s, score in enumerate(score_recognition_list):
        if score not in  df: continue
        if not same_axe:
            ax = fig.add_subplot(2,3,sp)
            sp += 1
            ax.set_title(score)
        
        l = []
        for i,index in enumerate(df.index):
            if type(index) == tuple:
                l.append(('{} '*len(index)).format(*index))
            else:
                l.append(index)
            #~ bottom = 0
            
            y = df[score+suffix][index]
            if same_axe:
                width =.5/6
                x = i+s*width
            else:
                x = i
                width =.8
                
            ax.bar(x, y, width =width, color = score_colors[score])
            #~ if y>0:
                #~ ax.text(i+.4, bottom+y/2, s,horizontalalignment='center',verticalalignment='center')
            #~ bottom += y
        if same_axe:
            ax.set_xticks(arange(len(l))+.25)
        else:
            ax.set_xticks(arange(len(l))+.4)
        ax.set_xticklabels(l, rotation = 45)


    
if __name__ =='__main__':
    #~ test_create_df_recognition_subject()
    #export_xl_recognition_scores_individual()
    
    #~ test_get_all()
    export_xl_recognition_scores_synthesis()
    #export_xl_timing_recognition_synthesis()
