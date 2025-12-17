"""
Script to compute Time-Frequency Power Features.

This script loads bipolarized data (`_bipo.npz`), computes power in specified
frequency bands (Delta to Gamma), and saves the results as `_power.npz` files
in the `feature/` folder.

Usage:
    Run directly: python 04_ComputePowerFeatures.py
    Ensures `config.py` is available.

Input:
    - `_bipo.npz` files (Bipolarized sEEG data)

Output:
    - `_power.npz` files containing:
        - `xpow`: Power matrix (n_elec, n_freqs, n_time_windows, n_trials)
        - `xpow_pval`: p-values (if computed)
        - `fname`: Frequency band names

Dependencies:
    - brainpipe.feature.power
"""
import numpy as np
from os import path
from time import sleep
import sys
from pathlib import Path

# Add project root to sys.path to find config.py
# olfaction_scripts_old is 2 levels up from 0_scripts_Etienne
current_dir = Path(__file__).resolve().parent
root_dir = current_dir.parent.parent
sys.path.append(str(root_dir))

try:
    import config
    STUDY_NAME = config.STUDY_NAME
except ImportError:
    STUDY_NAME = 'Olfacto'

from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt

if __name__ == "__main__":
    # Get all files include in the database
    try:
        st = study(STUDY_NAME)
    except NameError:
        print("Warning: 'study' not defined. Ensure brainpipe is installed.")
        sys.exit(1)

    # Search pattern needs to be generic or configurable? 
    # The original path 'database/TS_E_all_cond_by_block_trigs_filter1_500art/' seems very specific.
    # Refactoring to search in 'database/' globally or assume standard structure?
    # Keeping original folder logic but making it safer?
    # No, 'TS_E_all_cond_by_block_trigs_filter1_500art' implies a specific experiment run.
    # I will keep the folder string but acknowledge it in comments as 'Likely needs update'.
    
    target_folder = 'database/TS_E_all_cond_by_block_trigs_filter1_500art/'
    print(f"-> Searching for _bipo.npz in {target_folder}")
    files = st.search('_bipo.npz', folder=target_folder)
    
    if not files:
        print(f"-> No files found in {target_folder}. Checking root 'database/'...")
        files = st.search('_bipo.npz', folder='database/')
    
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
