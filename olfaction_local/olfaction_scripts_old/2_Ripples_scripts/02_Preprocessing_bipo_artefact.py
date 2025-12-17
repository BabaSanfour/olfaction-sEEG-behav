"""Script to bipolarize data.
"""
import numpy as np
from os import path

import pandas as pd
from mne.filter import resample, notch_filter, filter_data
from brainpipe.preprocessing import bipolarization
from brainpipe.system import study
from params import subjects
import matplotlib.pyplot as plt

from utils_functions import rename_elecs_df

def preprocessing_data(loadname):
    """
    Bipolarize all data using adjacent electrode
    Resample each bipolar derivation at 500 Hz
    Remove the 50-Hz power line interference (and its harmonics) using a zero-phase notch filter

    Parameters
    ----------
    file : npz files to concatenate
        Contains data ('x') and all electrodes information ('labels', 'channels', 'xyz')

    Returns
    -------
    bipo : an npz file with all data bipolarized
        Data of shape (nelecs x n_pts x ntrials).
    """
    mat = np.load(loadname, allow_pickle=True)

    # Bipolarize your data :
    x_bip, chan_bip, label_bip, xyz_bip = bipolarization(data=mat['x'],
                                            channel=[mat['channel'][i][0] for i in range(len(mat['channel']))],
                                            label=[mat['label'][i][0] for i in range(len(mat['label']))],
                                            xyz=mat['xyz'])

    # Filter and resample data
    x_bip = x_bip.swapaxes(1,2) #npts as the last dimension instead of trials
    x_bip = resample(x_bip,down=(512/500))
    x_bip = np.array(x_bip, dtype='float64')
    x_bip = notch_filter(x_bip, Fs=500., freqs=np.arange(50,250,50), method='fir')
    x_bip = x_bip.swapaxes(1,2)

    #rename elecs and remove elecs outside the brain
    columns = ['channels','labels','x','y','z']
    data = np.concatenate((np.array(chan_bip)[:,np.newaxis],
                            np.array(label_bip)[:,np.newaxis], xyz_bip),axis=-1)
    df = pd.DataFrame(data=data, columns=columns)
    new_df = rename_elecs_df(df)
    new_labels, new_chans = new_df['new_labels'].values, new_df['channels'].values
    print(np.unique(new_labels))
    idx_sel = [i for i,ch in enumerate(chan_bip) if ch in new_chans]

    # Save bipolarized and filtered data :
    mat = dict(mat)
    mat['x'], mat['channel'] = np.array(x_bip)[idx_sel], new_chans
    mat['xyz'], mat['label'] = xyz_bip[idx_sel], np.array(label_bip)[idx_sel]
    mat['new_label'] = np.array(new_labels)
    savename = loadname.replace('.npz', '_bipo.npz')
    np.savez(savename, **mat)

def detect_artefacts(loadnames):
    """
    Compute mean and sd on all electrodes and trials in order to detect potential artefacts in the data

    Parameters
    ----------
    files : list of filenames (monopolar or bipolarized signals) Encoding and Retrieval

    Returns
    -------
    detect_list : array
        Boolean array (nelecs, npts, ntrials) 1 representing potential artefacts
    """

    #take only common electrodes from encoding and retrieval in SELECTED ROIS
    mat_files = [np.load(loadname,allow_pickle=True) for loadname in loadnames]
    print([mat['x'].shape for mat in mat_files])
    print(np.unique(mat_files[0]['new_label']))
    labels_take = ['orbital','orbital2','HC','olf']
    list_lab = [[i for i,lab in enumerate(mat['new_label']) if lab in labels_take] for mat in mat_files]
    print('len ch sel', [len(x) for x in list_lab])
    channels = [mat['channel'][sel] for mat,sel in zip(mat_files,list_lab)]
    channels_common = list(set(channels[0]) & set(channels[1]) & set(channels[2]))
    print('nb common elecs', len(channels_common))

    #create index to identify common chans and trials
    idx_sels = [[i for i,chan in enumerate(chans) if chan in channels_common] for chans in channels]
    cum_trials = np.cumsum([mat['x'].shape[-1] for mat in mat_files])
    cum_trials = np.insert(cum_trials,0,0) #add a 0 to initiate indexes
    idx_trials = [[np.arange(cum_trials[i],cum_trials[i+1])] for i in range(len(cum_trials)-1)]

    #different nb of TIME POINTS (Rec longer so averaged for art detection)
    mean = np.concatenate(([np.mean(mat['x'][idx_sels[i]],axis=1) for i,mat in enumerate(mat_files)]),axis=-1)
    nchans, ntrials = mean.shape
    print('ALL data concat', mean.shape)
    chan_sel = [chans[idx_sels[i]] for i, chans in enumerate(channels)]
    #check if data is sorted in all arrays
    for a,b,c in zip(chan_sel[0],chan_sel[1],chan_sel[2]):
        if not a==b==c:
            print('warning !!!!! DATA not in the same ORDER')

    #detect artifact threshold according to Vaz et al. 2019
    Q1, Q3 = np.quantile(mean,0.25), np.quantile(mean,0.75)
    thr_e_tr = Q3+2.3*(Q3-Q1) #Vaz et al. 2019
    clean_id = mean < thr_e_tr #nchans
    print(np.sum(clean_id,axis=0),np.sum(clean_id,axis=1))

    #first remove chans with more than 5% of trials with artefacts
    thr_ch = ntrials * 0.98 if nchans > 2 else ntrials * 0.95
    good_chs = np.sum(clean_id,axis=1) > thr_ch
    clean_ch = clean_id[good_chs,:]
    nchans_good = clean_ch.shape[0]
    print('nb of chans removed', nchans - nchans_good)

    #second remove bad trials remaning in correct chans
    thr_tr = nchans_good
    good_tr = np.sum(clean_ch,axis=0) == thr_tr
    clean_tr = clean_ch[:,good_tr]
    print('nb of trials removed', ntrials - clean_tr.shape[-1])

    #select COMMON CHANS (idx_sel) with NO ARTEFACTS (good_ch)
    dict_ = {}
    dict_['new_label'] = mat_files[0]['new_label'][list_lab[0]][idx_sels[0]][good_chs]
    dict_['channel'] = chan_sel[0][good_chs]
    dict_['label'] = mat_files[0]['label'][list_lab[0]][idx_sels[0]][good_chs]
    print('new labels', np.unique(mat_files[0]['new_label'][list_lab[0]][idx_sels[0]][good_chs]))
    dict_['xyz'] = mat_files[0]['xyz'][list_lab[0]][idx_sels[0]][good_chs]

    #select good trials by exp condition
    for l,loadname in enumerate(loadnames):
        ind_sel = np.array(idx_sels[l],dtype=int)[good_chs]
        ind_tr = good_tr[tuple(idx_trials[l])]
        dict_['trials_sel'] = ind_tr
        for fi in ['odor','EM_2gr','EM_3gr']:
            dict_[fi] = mat_files[l][fi][ind_tr]
        data_sel = mat_files[l]['x'][ind_sel][:,:,ind_tr]
        print(reps[l], 'nb of removed trial', len(ind_tr) - np.sum(ind_tr))
        dict_['x'] = data_sel
        print('data to save',data_sel.shape)
        np.savez(loadname.replace('.npz','_clean.npz'), **dict_)

if __name__ == "__main__":

    process = 'artefact' #'preprocessing', 'artefact'
    reps =  ['E_odors', 'R_odors','R_rec']
    st = study('Ripples')

    if process == 'preprocessing':
        #Preprocess all files in the database (ignore files containing '_bipo'):
        for rep in reps:
            files = [k for k in st.search('_odor_', folder=('database/'+rep)) if not k.find('_bipo')+1]
            for fi in files:
                # Load file :
                loadname = path.join(st.path, 'database/'+rep, fi)
                print('-> Preprocess: '+loadname)
                preprocessing_data(loadname)

    if process == 'artefact':
        #Detect artefacts on bipo files
        for su in subjects:
            files = []
            for rep in reps:
                files.extend([k for k in st.search(su+'_cond=ALL_odors=ALL', folder=('database/'+rep)) if not k.find('_clean')+1])
            loadnames = [path.join(st.path, 'database/'+reps[i],fi) for i,fi in enumerate(files)]
            print('processing ',su)
            detect_artefacts(loadnames)
