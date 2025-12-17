"""This script will compute features on a variety
of frequency bands. This can be long, so may be you should
run it on a big computer. From the brainpipe folder feature,
I load three types of features :
    - amplitude (amplitude of hilbert transform)
    - power (can be either with wavelet / hilbert)
    - sigfilt (just the filtered signal)
"""
import numpy as np
from os import path
from time import sleep
from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt

if __name__ == "__main__":
    # Get all files include in the database
    st = study('Olfacto')
    files = st.search('_bipo.npz', folder='database/TS_E_all_cond_by_block_trigs_filter1_500art/')
    for fi in files:
        # Define power settings :
        kwargs = {} # Define an empty dictionnary to save all power parameters
        kwargs['f'] = [[2, 4], [4, 8], [8, 13], [13, 30], [30, 60], [60, 120]] # Frequency vector
        fname = ['delta', 'theta', 'alpha', 'beta', 'gamma30-60', 'gamma60-120'] # Name of each frequency
        #kwargs['baseline'] = (10, 1536) # Where your baseline (start, end) (IN SAMPLE) rest period ~ 3s
        #kwargs['norm'] = 0 # Type of normalisation (see help on power for more details)
        kwargs['width'], kwargs['step'] = 100, 50 # take power in 100 samples windows width every 50 samples
        # Load file :
        loadname = path.join(st.path, 'database/TS_E_all_cond_by_block_trigs_filter1_500art/', fi)
        print('-> Compute power on: '+loadname)
        mat = np.load(loadname)
        x, sf = mat['x'], int(mat['sf'])
        n_elec, n_pts, n_trials = x.shape
        # Define a power object :
        powObj = power(sf, n_pts, **kwargs)
        kwargs['xpow'],  kwargs['xpow_pval']= powObj.get(x, n_jobs=-1)
        # Finally save it :
        kwargs['fname'] = fname
        savename = loadname.replace('.npz', '_power.npz').replace('database', 'feature')
        np.savez(savename, **kwargs)
        del kwargs, powObj, mat, x, sf, n_elec, n_trials, n_pts, fname
