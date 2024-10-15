# -*- coding: utf-8 -*-

from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

from respiration_features import create_df_respiration_features_subject
from score_episodic import create_df_episodic_subject, add_tr_by_score

from collections import OrderedDict

from matplotlib import pyplot


l = ['WWW', 'WWhich', 'WWhere', 'What' ]
score_episodic_strict_list =  [ ]
for e in l:
    score_episodic_strict_list.append(e+'-S1')


def get_all():
    all = [ ]  
    for subject in session.query(Subject).order_by('Subject.`group`'):
        df1 = create_df_episodic_subject(subject)
        df2 = create_df_respiration_features_subject(subject, 'R', 'total')
        df_1_2= concat([df1, df2], axis = 1)
        #~ print df_1_2.columns
        all.append(df_1_2)
    #~ print all
    return concat(all, axis=0)
    
def add_respi_by_score(all, select, df):
    prefix = 'total_'
    resp_features = ['mean_insp_duration', 'mean_insp_volume', 'mean_insp_max',
                        'mean_exp_duration', 'mean_exp_volume', 'mean_exp_max',
                        'mean_cycle_duration','nb_cycle', 'mean_frequency',
                        'inv_mean_duration', 
                        ]
    for e in resp_features:
        for s in score_episodic_strict_list:
            df[e+'_'+s] = Series(index= df.index)
    for name, subdf in select:
        for e in resp_features:
            for s in score_episodic_strict_list:
                df.loc[name, e+'_'+s] = subdf.loc[subdf['score_episodic_strict'] == s, prefix+e].mean()

def export_episodic_respi_synthesis():
    all = get_all()
    writer = ExcelWriter('Respi_epi_synthesis.xls')
    
    # par sujet strict
    select = all.groupby(['subject_name'])
    score_by_subject1= select[score_episodic_strict_list].sum() 
    add_tr_by_score (all, select, score_by_subject1)
    add_respi_by_score(all, select, score_by_subject1)
    score_by_subject1.to_excel(writer, sheet_name='by_subject_strict')
    print score_by_subject1

    # par context_day strict
    select = all.groupby(['subject_name', 'context_day'])
    score_by_context_day1= select[score_episodic_strict_list].sum() 
    add_tr_by_score (all, select, score_by_context_day1)
    add_respi_by_score(all, select, score_by_context_day1)
    score_by_context_day1.to_excel(writer, sheet_name='by_context_day_strict')
    print score_by_context_day1
    
    # par image context strict
    select = all.groupby(['subject_name', 'image_context'])
    score_by_image_context1= select[score_episodic_strict_list].sum() 
    add_tr_by_score (all, select, score_by_image_context1)
    add_respi_by_score(all, select, score_by_image_context1)
    score_by_image_context1.to_excel(writer, sheet_name='by_image_context_strict')
    print score_by_image_context1

    
    all.to_excel(writer, sheet_name='all')
    
    writer.save()
    
if __name__ =='__main__':

    export_episodic_respi_synthesis()
    #~ get_all()