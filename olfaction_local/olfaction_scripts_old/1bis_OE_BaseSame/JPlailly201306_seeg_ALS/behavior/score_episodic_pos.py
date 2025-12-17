# -*- coding: utf-8 -*-

from connection import *
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

from respiration_tools import get_time_first_inspi_or_trig
import codes
from score_common import create_df_subject
from score_reconnaissance import create_df_recognition_subject

from collections import OrderedDict
from matplotlib import pyplot

pos_files_dirname = './FichiersPOS'
if not os.path.exists(pos_files_dirname):
    os.mkdir(pos_files_dirname)


def codes_triggers(trial):
        
    score1 = None
       
    if trial.score_recognition == 'miss':
        score1 = 2
    if trial.score_recognition == 'cr':
        score1 = 3
    if trial.score_recognition == 'fa':
        score1 = 4
    if trial.score_recognition == 'hit':
        image_pos, spot_pos,spot_pos_group,  image = codes.num_odor_to_ref[trial.run.subject.name][trial.odor_num]
        if trial.odor_num not in codes.num_odor_to_ref[trial.run.subject.name]:
            return None, None
        elif (trial.selected_image ==image_pos) and (trial.click_posistion==spot_pos):
            #WWW
            score1 = 15
        elif (trial.selected_image ==image_pos) and (trial.click_posistion!=spot_pos):
            #WWhich
            score1 = 16
        elif (trial.selected_image !=image_pos) and (trial.click_posistion==spot_pos):
            #WWhere
            score1 = 17
        elif (trial.selected_image !=image_pos) and (trial.click_posistion!=spot_pos):
            #What
            score1 = 18
        else:
            print 'BUGGGGGGG'
    
    return score1
    
def create_df_retrieval_pos(subject):
    ix = [ ]
    for run in subject.runs.filter_by(exp='R').order_by('run_index'):
        for trial in run.trials:
            ix.append(trial.id)
    n = len(ix)
    df = DataFrame(index = ix)
    
    df['retrieval_session']  = Series(index = ix, dtype = 'S')
    for k in ['start', 'odor', 'first_inspi', 'rec', 'epi']:
        df[k] = Series(index = ix, dtype = int)
        df[k+'_code'] = Series(index = ix, dtype = int)

    for run in subject.runs.filter_by(exp = 'R'):
        resp = run.respirationsignals[0]    
        for trial in run.trials:
            score1 = codes_triggers(trial)
            #~ if score1 is not None:
            df.loc[trial.id,'retrieval_session'] = run.index
            df.loc[trial.id,'start']= int(trial.time*512.)
            df.loc[trial.id,'start_code'] = score1+100

            df.loc[trial.id,'odor'] = int(trial.triggered_odor_time*512.)
            df.loc[trial.id,'odor_code'] = score1+200
            
            if trial.first_inspiration_time != 0.:
                df.loc[trial.id,'first_inspi'] = int(trial.first_inspiration_time*512.)
            else: df.loc[trial.id,'first_inspi'] = int(trial.triggered_odor_time*512.)
            df.loc[trial.id,'first_inspi_code'] = score1+300
            
            df.loc[trial.id,'rec'] = int((trial.time+trial.recognition_time)*512.)
            df.loc[trial.id,'rec_code'] = score1+400
            
            if trial.context_time != 0.:
                df.loc[trial.id,'epi']= int((trial.time+trial.context_time)*512.)
                df.loc[trial.id,'epi_code'] = score1+500
            else: continue
    session.commit()
    
    return df

def export_retrieval_pos_all():
    
    writer = ExcelWriter('Fichiers_POS_Retrieval.xls')
    for subject in session.query(Subject).order_by('Subject.`group`'):
        print subject.name
        df =  create_df_retrieval_pos(subject)
        df.to_excel(writer, sheet_name=subject.name)
        
        
        for retrieval_session in [1,2,3]:
            all = []
            for k in ['start', 'odor', 'first_inspi', 'rec', 'epi']:
                df2 = df[df['retrieval_session']==retrieval_session]
                df2 = df2[[k, k+'_code']]
                df2.columns = ['sample', 'code']
                df2['reject'] = 0
                all.append(df2)
            all = concat(all, axis=0)
            all = all[~all['sample'].isnull()]
            all.to_csv(os.path.join(pos_files_dirname,'{}_{}{}.pos'.format(subject.name, 'R', retrieval_session)), header = False, index = False, sep = '\t', float_format = '%.0f' )
        
        
        
        
        
    writer.save()
    
    
    

   
if __name__ =='__main__':
    export_retrieval_pos_all()
    
