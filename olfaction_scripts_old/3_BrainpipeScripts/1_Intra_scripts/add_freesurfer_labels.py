import numpy as np
import pandas as pd
from utils import subjects
from os.path import join
from brainpipe.system import study

st = study('Olfacto')
FSF_file = join(st.path, 'feature_new/All_subjects_Freesurfer_labels.csv')

def add_fsf_labels(df):
    df_fsf = pd.read_csv(FSF_file)
    free_labs, hip_CA, hip_fs60 = [], [], []
    for su in subjects:
        df_su = df.loc[df['subjects']==su]
        if df_su.shape[0]:
            df_f_su = df_fsf.loc[df_fsf['subjects']==su]
            list_chans = [lab.replace("'","") for lab in df_su['channels'].str.upper()] \
                                        if su == 'CHAF' else df_su['channels'].str.upper().values
            idx = [i for i,lab in enumerate(df_f_su['contact'].values) if lab in list_chans]

            #due to contact type 3x5 some are missing in Freesurfer data
            missing = [elec for elec in list_chans if elec not in df_f_su['contact'].values[idx]]
            idx_miss = [i for i,elec in enumerate(list_chans) if elec in missing]
            idx_miss = [ind -i for i, ind in enumerate(idx_miss)]
            print(su, 'missing contacts', missing, idx_miss)

            #select all subject's labels
            fr_l = df_f_su['Freesurfer'].values[idx]
            ca_l = df_f_su['hip_CA'].values[idx]
            fs_l = df_f_su['hip_FS60'].values[idx]

            #add None labels at the good index location in subject's df
            fr_l = np.insert(fr_l, idx_miss, 'not_fs')
            ca_l = np.insert(ca_l, idx_miss, 'not_fs')
            fs_l = np.insert(fs_l, idx_miss, 'not_fs')
            free_labs.extend(fr_l), hip_CA.extend(ca_l), hip_fs60.extend(fs_l)

    df['fsf'], df['hip_CA'], df['hip_FS60'] = free_labs, hip_CA, hip_fs60
    return df
