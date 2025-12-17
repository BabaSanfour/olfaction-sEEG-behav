# -*- coding: utf-8 -*-

from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex
#~ from respiration_features import get_respiration_feature


from score_common import create_df_subject
from respiration_tools import get_time_first_inspi_or_trig
from respiration_features import get_respiration_feature_for_trial
import codes

from collections import OrderedDict

from matplotlib import pyplot


l = ['WWW', 'WWhich', 'WWhere', 'What' ]


score_episodic_strict_list =  [ ]
score_episodic_elargi_list = [ ]
for e in l:
    score_episodic_strict_list.append(e+'-S1')
    score_episodic_elargi_list.append(e+'-S4')

score_colors  = { }
for  e in score_episodic_strict_list+score_episodic_elargi_list:
    e1,e2 = e.split('-')
    colors =  {
                        'WWW' : '#9A2EFE',
                        'WWhich' : '#BCA9F5',
                        'WWhere' : '#BCA9F5', 
                        'What'  :  '#BCA9F5', 
                        }
    score_colors[e] = colors[e1]
    #~ print e, prefix



filename_results1 = 'episodic_individual_results_.xls'
filename_results2 = 'episodic_summary.xls'
filename_results3 = 'episodic_ind_by_rep.xls'



#~ show_figures = True
show_figures = False


def score_episodic(trial):
    
    if trial.odor_num not in codes.num_odor_to_ref[trial.run.subject.name]:
        return None, None
    
    score1, score2 = None, None
    image_pos, spot_pos,spot_pos_group,  image = codes.num_odor_to_ref[trial.run.subject.name][trial.odor_num]
    
    if trial.score_recognition == 'hit':
        if (trial.selected_image ==image_pos) and (trial.click_posistion==spot_pos):
            score1 = 'WWW-S1'
        elif (trial.selected_image ==image_pos) and (trial.click_posistion!=spot_pos):
            score1 = 'WWhich-S1'
        elif (trial.selected_image !=image_pos) and (trial.click_posistion==spot_pos):
            score1 = 'WWhere-S1'
        elif (trial.selected_image !=image_pos) and (trial.click_posistion!=spot_pos):
            score1 = 'What-S1'
        else:
            print 'BUGGGGGGG'
    
    if trial.score_recognition == 'hit':
        if (trial.selected_image ==image_pos) and (trial.click_posistion in spot_pos_group):
            score2 = 'WWW-S4'
        elif (trial.selected_image ==image_pos) and (trial.click_posistion not in spot_pos_group):
            score2 = 'WWhich-S4'
        elif (trial.selected_image !=image_pos) and (trial.click_posistion in spot_pos_group):
            score2 = 'WWhere-S4'
        elif (trial.selected_image !=image_pos) and (trial.click_posistion not in spot_pos_group):
            score2 = 'What-S4'
        else:
            print 'BUGGGGGGG'
    
    if trial.score_recognition == 'miss':
        score1, score2 = 'miss', 'miss'
    
    return score1, score2
    
    
def create_df_episodic_subject(subject):
    df =  create_df_subject(subject)
    ix = df.index
    ## ajjout de colonne
    df['score_episodic_strict']  = Series(index = ix, dtype = 'S')
    df['score_episodic_elargi']  = Series(index = ix, dtype = 'S')
    df['first_inspi_or_trig']  = Series(index = ix)
    df['odor_recognition_rt']  = Series(index = ix)
    df['episodic_rt']  = Series(index = ix)
    for s in score_episodic_strict_list:
        df[s] = Series(index = ix)
    for s in score_episodic_elargi_list:
        df[s] = Series(index = ix)
    
    
    # remplissage
    for run in subject.runs.filter_by(exp = 'R'):
        if run.trials.count()==0:continue
        resp = run.respirationsignals[0]
        resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
        for trial in run.trials:
            
            score1, score2 = score_episodic(trial)
            
            trial.score_episodic_strict = df.loc[trial.id,'score_episodic_strict'] = score1
            trial.score_episodic_elargi = df.loc[trial.id,'score_episodic_elargi'] = score2
            
            trial.score_episodic_strict = score1
            trial.score_episodic_elargi = score2
            
            for s in score_episodic_strict_list:
                df.loc[trial.id,s] = int(s == score1)            
            for s in score_episodic_elargi_list:
                df.loc[trial.id,s] = int(s == score2)     
        
            df.loc[trial.id,'first_inspi_or_trig'] = get_time_first_inspi_or_trig(trial, resp, max_percent = max_percent_for_respiration_cycle)
            if df.loc[trial.id,'first_inspi_or_trig'] is not None:
                df.loc[trial.id,'odor_recognition_rt'] = trial.recognition_time - (df.loc[trial.id,'first_inspi_or_trig'] -trial.time)
                df.loc[trial.id,'episodic_rt'] = trial.context_time - (df.loc[trial.id,'first_inspi_or_trig'] -trial.time)

            prefix = 'episodic_'
            t1 = trial.recognition_time+trial.time 
            t2 = trial.time+25.
            features = get_respiration_feature_for_trial(trial, resp_f, t1, t2, n_cycle = 15)
            for k, v in features.items():
                if prefix+k not in df.columns:
                    df[prefix+k] = Series(index = ix)
                df.loc[trial.id,prefix+k] = v
    subject.score_episodic_wwwS1 = float(sum(df['score_episodic_strict']=='WWW-S1'))
    #~ subject.score_episodic_wwwS4 = float(sum(df['score_episodic_elargi']=='WWW-S4'))
    session.commit()
    return df
    
def test_create_df_episodic_subject():
    subject = session.query(Subject).filter_by(name = 'LEFC').first()
    df =  create_df_episodic_subject(subject)
    print df
    print df[ ['is_target','image_context', 'score_episodic_strict', 'score_episodic_elargi', 'trial_click_posistion']]


def export_xl_episodic_scores_individual():
    
    writer = ExcelWriter(filename_results1)
    
    for subject in session.query(Subject).order_by('Subject.`group`'):
    #~ for subject in session.query(Subject).order_by('Subject.`group`').filter_by(name = 'LEFC'):
        print subject.name
        df = create_df_respiration_features_subject(subject, 'R', 'episodic')
        #df =  create_df_episodic_subject(subject)
        df.to_excel(writer, sheet_name=subject.name)
        #df2 = df[['subject_name', 'odor_num', 'score_recognition']]
        #df2.to_excel(writer, sheet_name=subject.name)
    writer.save()

def export_xl_episodic_ind_rep0():
    
    writer = ExcelWriter(filename_results3)
    
    for subject in session.query(Subject).order_by('Subject.`group`'):
    #~ for subject in session.query(Subject).order_by('Subject.`group`').filter_by(name = 'LEFC'):
        print subject.name
        df =  create_df_episodic_subject(subject)
        select = df.groupby(['odor_num'])
        df2 = select['subject_name','odor_num','repetition_of_target','score_episodic_strict','score_episodic_elargi']
        print df2
        df2.to_excel(writer,sheet_name=subject.name)
        #~ df2 = df[['subject_name', 'odor_num', 'score_recognition']]
        #~ df2.to_excel(writer, sheet_name=subject.name)
    writer.save() 


def get_all():
    all = [ ] 
    from respiration_features import create_df_respiration_features_subject
    for subject in session.query(Subject).order_by('Subject.`group`'):
        df = create_df_episodic_subject(subject)
        all.append(df)
    return concat(all)

def test_get_all():
    all = get_all()
    print all

def add_tr_by_score(all, select, df):
    for s in score_episodic_strict_list:
        df[s+'_odor_rt'] = Series(index = df.index)
        df[s+'_context_rt'] = Series(index = df.index)
        
    for name, group in select:
        for s in score_episodic_strict_list:
            sub = group[group['score_episodic_strict'] == s]
            df[s+'_odor_rt'][name] = sub['odor_recognition_rt'].mean()
            df[s+'_context_rt'][name] = sub['episodic_rt'].mean()

def export_xl_episodic_scores_synthesis():
    all = get_all()
    
    writer = ExcelWriter(filename_results2)
    
    # par sujet strict
    all2 = all[all['is_target'] ==1]
    select = all2.groupby(['subject_name'])
    #print all2.columns
    score_by_subject1= select[score_episodic_strict_list].sum()
    add_tr_by_score (all, select, score_by_subject1)
    score_by_subject1.to_excel(writer, sheet_name='by_subject')
    #print score_by_subject1

    # par sujet strict
    all2 = all[all['is_target'] ==1]
    select = all2.groupby(['subject_name','odor_num'])
    print select
    score_by_subject1= select[score_episodic_strict_list].sum()
    add_tr_by_score (all, select, score_by_subject1)
    score_by_subject1.to_excel(writer, sheet_name='by_odor')
    #print score_by_subject1

    # par context_day strict
    select = all.groupby(['subject_name', 'context_day'])
    score_by_context_day1= select[score_episodic_strict_list].sum() 
    add_tr_by_score (all, select, score_by_context_day1)
    score_by_context_day1.to_excel(writer, sheet_name='by_context_day_strict')
    #print score_by_context_day1
    
    # par image context strict
    select = all.groupby(['subject_name', 'image_context'])
    score_by_image_context1= select[score_episodic_strict_list].sum() 
    add_tr_by_score (all, select, score_by_image_context1)
    score_by_image_context1.to_excel(writer, sheet_name='by_image_context_strict')
    #print score_by_image_context1


    
    # par sujet elargi
    #~ score_by_subject2= all[['subject_name']+score_episodic_elargi_list].groupby(['subject_name']).sum()
    #~ score_by_subject2.to_excel(writer, sheet_name='by_subject_elargi')
    #~ print score_by_subject2
   
    # par context_day elargi
    #~ score_by_context_day2= all[['subject_name', 'context_day', ]+score_episodic_elargi_list].groupby(['subject_name', 'context_day']).sum()
    #~ score_by_context_day2.to_excel(writer, sheet_name='by_context_day_elargi')
    #~ print score_by_context_day2

    # par image context elargi
    #~ score_by_image_context2 = all[['subject_name', 'image_context', ]+score_episodic_elargi_list].groupby(['subject_name', 'image_context']).sum()
    #~ score_by_image_context2.to_excel(writer, sheet_name='by_image_context_elargi')
    #~ print score_by_image_context2
    
    
    writer.save()
    
    if show_figures:
        plot_df(score_by_subject1, score_episodic_strict_list, title = 'subject')
        #~ plot_df(score_by_subject2, score_episodic_elargi_list, title = 'subject')
        #~ plot_df(score_by_odor1, score_episodic_strict_list, title = 'subject/odor')
        #~ plot_df(score_by_odor2, score_episodic_elargi_list, title = 'subject/odor')
        #~ plot_df(score_by_repetition_of_target1, score_episodic_elargi_list)
        #~ plot_df(score_by_repetition_of_target2, score_episodic_elargi_list)
        #~ plot_df(score_by_context_day1, score_episodic_strict_list, title = 'subject/context_day')
        #~ plot_df(score_by_context_day2, score_episodic_elargi_list, title = 'subject/context_day')
        #~ plot_df(score_by_image_context1, score_episodic_strict_list, title = 'subject/image_context')
        #~ plot_df(score_by_image_context2, score_episodic_elargi_list, title = 'subject/image_context')
        
        pyplot.show()



def plot_df(df,score_list, title = ''):

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
        for s in score_list:
            if s in df:
                y = df[s][index]
                ax.bar(i, y, width =.8, color = score_colors[s], bottom = bottom)
                if y>0:
                    ax.text(i+.4, bottom+y/2, s,horizontalalignment='center',verticalalignment='center')
                bottom += y
    ax.set_xticks(arange(len(l))+.4)
    ax.set_xticklabels(l)
    

    
if __name__ =='__main__':
    #~ test_create_df_episodic_subject()
    #export_xl_episodic_scores_individual()
    #export_xl_episodic_ind_rep0()
    
    #~ test_get_all()
    export_xl_episodic_scores_synthesis()
