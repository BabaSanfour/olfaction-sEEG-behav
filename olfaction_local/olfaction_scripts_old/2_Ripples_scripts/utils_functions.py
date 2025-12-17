import numpy as np
from os import listdir
from scipy.signal import convolve2d

###################################################
x0,y0,z0 = 20, 30,-16 #VOI for olf OFC
x1,y1,z1 = -20, 32,-16 #VOI for olf OFC
x0_ins,y0_ins,z0_ins = 42.0, -7.3, 0.3 #VOI for olf Ins
x1_ins,y1_ins,z1_ins = -38, -4, 12 #VOI for olf Ins
rad = 11 #radius of the VOI
HC_lim = -20.5 #limit anterior/posterior part of the aHC
###################################################

def revert_dico(dico):
    new_dico = {}
    for su in dico:
        dico_su = {}
        for cond in dico[su]:
            dico_su.update({od:cond for od in dico[su][cond]})
        new_dico[su] = dico_su
    return new_dico

def conv2(x, y, mode='same'):
    """Equivalent of conv2 MATLAB_R2019 function"""
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)

def rename_elecs(names,x,y,z):

    labels = ['pPirT' if lab =='Amg-PirT' else lab for lab in names]

    id_R_OFC = [1 if all([x0-rad<=round(x[i])<=x0+rad, y0-rad<=round(y[i])<=y0+rad,
                    z0-rad<=round(z[i])<=z0+rad]) else 0 for i in range(len(labels))]
    id_L_OFC = [1 if all([x1-rad<=round(x[i])<=x1+rad, y1-rad<=round(y[i])<=y1+rad,
                    z1-rad<=round(z[i])<=z1+rad]) else 0 for i in range(len(labels))]
    idx_OFC = np.array(id_R_OFC) + np.array(id_L_OFC)
    labels = ['OFC_olf' if i == 1 else lab for i,lab in zip(idx_OFC,labels)]

    id_R_ins = [1 if all([x0_ins-rad<=round(x[i])<=x0_ins+rad, y0_ins-rad<=round(y[i])<=y0_ins+rad,
                    z0_ins-rad<=round(z[i])<=z0_ins+rad]) else 0 for i in range(len(labels))]
    id_L_ins = [1 if all([x1_ins-rad<=round(x[i])<=x1_ins+rad, y1_ins-rad<=round(y[i])<=y1_ins+rad,
                    z1_ins-rad<=round(z[i])<=z1_ins+rad]) else 0 for i in range(len(labels))]
    idx_ins = np.array(id_R_ins) + np.array(id_L_ins)
    labels = ['Ins_olf' if i == 1 else lab for i,lab in zip(idx_ins,labels)]
    return np.array(labels)

def rename_elecs_df(df):
    new_dtypes = {"x": np.float64, "y": np.float64, "z":np.float64}
    df = df.astype(new_dtypes)
    names = rename_elecs(df['labels'],df['x'],df['y'],df['z'])
    names_up = [x.upper() for x in names] #homogenize names
    y = df['y'].values

    matches_parietal = ['SMG','ANG','PCUN','PRG','POG','RSC']
    matches_MTL = ['AHC','MHC','PHC']
    new_labels = ['orbital' if 'ORG' in lab else lab for lab in names_up]
    new_labels = ['orbital' if 'OFC_OLF' in lab else lab for lab in new_labels]
    new_labels = ['olf' if any(x in lab for x in ['AMG','PIR','ENT']) else lab for lab in new_labels]
    new_labels = ['HC' if any(x in lab for x in matches_MTL) else lab for lab in new_labels]
    new_labels = ['PHG_FuG' if any(x in lab for x in ['FUG','PHG']) else lab for lab in new_labels]
    new_labels = ['orbital2' if any(x in lab for x in ['IFG','ROG','IFPG']) else lab for lab in new_labels]
    new_labels = ['INS' if any(x in lab for x in ['INS','IG']) else lab for lab in new_labels]
    new_labels = ['temp_cx' if any(x in lab for x in ['TG','TP','PPO']) else lab for lab in new_labels]
    new_labels = ['cingulate' if 'CC' in lab else lab for lab in new_labels]
    new_labels = ['Frontal' if 'FG' in lab else lab for lab in new_labels]
    new_labels = ['Parietal' if any(x in lab for x in matches_parietal) else lab for lab in new_labels]
    new_labels = ['occipital' if any(x in lab for x in ['OCG','17']) else lab for lab in new_labels]
    new_labels = ['aMTL' if lab == 'MTL' and y[i] > -26 else lab \
                                  for i,lab in enumerate(new_labels)]
    df['new_labels'] = new_labels

    df2 = df[~df.new_labels.str.contains("SKULL")] #clean eletrodes name
    df2 = df2[~df2.new_labels.str.contains("LCR")] #clean eletrodes name
    df2 = df2[~df2.new_labels.str.contains("WM")] #clean eletrodes name
    df2 = df2[~df2.new_labels.str.contains("OUT")] #clean eletrodes name
    return df2

def compute_big_npz(rep, st_search,val_concat,val_unique):
    """
    Compute big files with all subjects' data (averaged across trials)

    Parameters
    ----------
    rep : string (path to the data)
    st_search : string (string to look for in the repository)
    val_concat : string (arrays to concatenate across patients)
    val_unique : string (arrays to include in the file - information about the analysis)

    Return
    ------
    big_npz : dict containing all features and electrodes information TO BE SAVED
    """

    big_npz = {}
    files = [f for f in listdir(rep) if (f.endswith(st_search)) & ('All_subj' not in f)]
    big_npz['subjects'] = np.concatenate([[fi.split('_')[0]]*np.load(rep+fi)['channel'].shape[0] for fi in files])

    for v in val_unique:
        big_npz[v] = np.load(rep+files[0],allow_pickle=True)[v]

    for v_c in val_concat:
        concat_dat = np.array([])
        for f in files:
            data = np.load(rep+f,allow_pickle=True)[v_c]
            if len(data.shape) < 3:
                concat_dat = np.concatenate((concat_dat,data),axis=0) if np.size(concat_dat) else data
            elif len(data.shape) > 3:
                data_mean = np.mean(data, axis=-1) #average across trials
                concat_dat = np.concatenate((concat_dat,data_mean),axis=1) if np.size(concat_dat) else data_mean
        big_npz[v_c] = concat_dat
    return big_npz
