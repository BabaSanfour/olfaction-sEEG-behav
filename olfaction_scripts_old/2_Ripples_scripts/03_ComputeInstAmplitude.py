"""This script will compute features on a variety of frequency bands.
From the brainpipe folder feature, I load three types of features :
    - amplitude (amplitude of hilbert transform)
    - power (can be either with wavelet / hilbert)
    - sigfilt (just the filtered signal)
"""

import numpy as np
from os.path import join, exists
from os import makedirs
from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt
from brainpipe.feat.utils._feat import _manageWindow
from params import subjects
from itertools import product

def compute_features(mat,x, norm=False):
    """
    Compute signal features (amplitude, power and sig filt)
    Compute statistics with baseline (data normalized or not)
    Save all data in an npz file
    Ripples borns and HFA (Vaz et al. 2019, 2020)

    Parameters
    ----------
    mat : npz files to process
        Contains data ('x') and all electrodes information ('labels', 'channels', 'xyz')
    x : array (nelecs, npts, trials)

    Return
    ------
    kwargs : dict containing all features and electrodes information TO BE SAVED
    """

    # Define features settings :
    kwargs = {} # Define an empty dictionnary to save all power parameters
    kwargs['f'] = [[80,120],[15,60],[70,200]] # Frequency vector
    fname = ['ripple', 'low_freq','HFA'] # Name of each frequency
    kwargs['baseline'] = None # (start, end) (IN SAMPLE) rest period ~ 3s (for stats: 1s before onset)
    kwargs['split'] = [4,None,13]
    kwargs['norm'] = 3 if norm == True else 0 # Type of normalisation (see help on power for more details)
    sf = 500
    kwargs['width'], kwargs['step'] = None, None #instantaneous amplitude

    # Load file :
    n_elec, n_pts, n_trials = x.shape

    # Define and compute feature objects :
    amplObj = amplitude(sf, n_pts, **kwargs)
    _, time = _manageWindow(n_pts, width=None, step=None)
    kwargs['time'] = np.array(time) / sf
    kwargs['xampl'], kwargs['pval_xampl'] = amplObj.get(x, statmeth='wilcoxon',n_jobs=-1)

    # Finally save it :
    kwargs['x'] = data
    kwargs['fname'], kwargs['new_labels'] = fname, mat['new_labels']
    kwargs['channels'], kwargs['labels'], kwargs['xyz'] = mat['channels'], mat['labels'], mat['xyz']
    return kwargs


if __name__ == "__main__":

    st = study('Ripples')
    reps = ['Retrieval_new_odors/', 'Retrieval_new_rec/','Encoding/']
    steps = ['']

    for rep, step in product(reps,steps):
        files = st.search('INFO.npz', folder=('database/'+rep))
        rep_save = join(st.path, 'feature/'+rep)
        if not exists(rep_save):
            makedirs(rep_save)
        for fi in files:
            print('Processing --> ', rep, fi, step)
            loadname = join(st.path, 'database/'+rep+fi)
            mat = np.load(loadname, allow_pickle=True)
            data = mat['x']
            if step == '':
                dico = compute_features(mat,data,norm=False)
                savename = loadname.replace('.npz', '_ampl.npz').replace('database', 'feature')
            elif step == 'Norm':
                dico = compute_features(mat,data,norm=True)
                savename = loadname.replace('.npz', '_feat_norm.npz').replace('database', 'feature')
            np.savez(savename, **dico)
