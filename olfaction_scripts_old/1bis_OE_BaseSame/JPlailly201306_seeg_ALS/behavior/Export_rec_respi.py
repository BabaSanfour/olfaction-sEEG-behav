# -*- coding: utf-8 -*-

from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

from respiration_features import create_df_respiration_features_subject
from score_reconnaissance import create_df_recognition_subject, add_tr_by_score

from collections import OrderedDict

from matplotlib import pyplot



score_recognition_list = [ 'hit',  'miss', 'cr', 'fa' ]

def get_all():
    all = [ ]  
    for subject in session.query(Subject).order_by('Subject.`group`'):
        df1 = create_df_recognition_subject(subject)
        df2 = create_df_respiration_features_subject(subject, 'R', 'recognition')
        df_1_2= concat([df1, df2], axis = 1)
        #~ print df_1_2.columns
        all.append(df_1_2)
    #~ print all
    return concat(all, axis=0)
    
def add_respi_by_score(all, select, df):
    prefix = 'recognition_'
    resp_features = ['mean_insp_duration', 'mean_insp_volume', 'mean_insp_max',
                        'mean_exp_duration', 'mean_exp_volume', 'mean_exp_max',
                        'mean_cycle_duration','nb_cycle', 'mean_frequency',
                        'inv_mean_duration', 
                        ]
    for e in resp_features:
        for s in score_recognition_list:
            df[e+'_'+s] = Series(index= df.index)
    for name, subdf in select:
        for e in resp_features:
            for s in score_recognition_list:
                df.loc[name, e+'_'+s] = subdf.loc[subdf['score_recognition'] == s, prefix+e].mean()

def export_recognition_respi_synthesis():
    all = get_all()
    writer = ExcelWriter('Respi_rec_synthesis.xls')
    
    # par sujet
    select = all.groupby(['subject_name'])
    score_by_subject = select[score_recognition_list].sum()
    add_tr_by_score(all, select, score_by_subject)
    add_respi_by_score(all, select, score_by_subject)
    print score_by_subject
    score_by_subject.to_excel(writer, sheet_name='by_subject')

    # par context_day
    all2 = all[all['is_target'] ==1]
    select = all2.groupby(['subject_name', 'context_day'])
    score_by_context_day = select['hit', 'miss'].sum()
    add_tr_by_score (all, select['hit', 'miss'], score_by_context_day)
    add_respi_by_score(all, select['hit', 'miss'], score_by_context_day)
    score_by_context_day.to_excel(writer, sheet_name='by_context_day')
    print  score_by_context_day
    
    # par image context
    all2 = all[all['is_target'] ==1]
    select = all2.groupby(['subject_name', 'image_context'])
    score_by_image_context = select['hit', 'miss'].sum()
    add_tr_by_score (all, select['hit', 'miss'], score_by_image_context)
    add_respi_by_score(all, select['hit', 'miss'], score_by_image_context)
    score_by_image_context.to_excel(writer, sheet_name='by_image_context')
    print score_by_image_context
    
    all.to_excel(writer, sheet_name='all')
    
    writer.save()
    
if __name__ =='__main__':

    export_recognition_respi_synthesis()
    #~ get_all()