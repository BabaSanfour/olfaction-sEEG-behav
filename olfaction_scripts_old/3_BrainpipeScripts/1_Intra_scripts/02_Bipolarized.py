"""
Script to bipolarize sEEG data.

This script loads preprocessed .npz files (cleaned data), applies bipolar referencing
(re-referencing each channel to its neighbor), and saves the result as `_bipo.npz`.

Usage:
    Run directly: python 02_Bipolarized.py
    Ensures `config.py` is available in the root.

Dependencies:
    - brainpipe.preprocessing.bipolarization
    - brainpipe.system.study
"""
import numpy as np
from os import path

from brainpipe.preprocessing import bipolarization
from brainpipe.system import study


import sys
from pathlib import Path

# Add project root to sys.path to find config.py
# olfaction_scripts_old is 2 levels up from 1_Intra_scripts
current_dir = Path(__file__).resolve().parent
root_dir = current_dir.parent.parent
sys.path.append(str(root_dir))

try:
    import config
    STUDY_NAME = config.STUDY_NAME
    # If using config, study path is handled by config + brainpipe logic
except ImportError:
    STUDY_NAME = 'Ripples' # Default fallback or 'Olfacto' as per other scripts? 
    # Note: 00_CreateStudy uses 'Olfacto', this file had 'Ripples'. 
    # Sticking to 'Ripples' for now if hardcoded, or config if available.
    # Actually, config.py defines STUDY_NAME = 'Olfacto'. 
    # If this script specifically needs 'Ripples', strictly overriding might break logic if study names matter.
    # However, 'Ripples' likely refers to a different project or copy-paste error? 
    # Given the directory is '1_Intra_scripts' and not '2_Ripples_scripts', 'Olfacto' is more likely correct.
    # Let's use config.STUDY_NAME if import succeeds.

if __name__ == "__main__":
    
    # Define subfolders to search in
    reps =  [   
            'Retrieval_new_odors',
            'Retrieval_new_rec',
        ]
    
    # Initialize Study
    # Note: If brainpipe is missing, this will fail.
    try:
        st = study(STUDY_NAME)
    except NameError:
        print("Warning: 'study' not defined. Ensure brainpipe is installed.")
        st = None
        
    if st:

for rep in reps:
    files = [k for k in st.search('_odor_', folder=('database/'+rep)) if not k.find('_bipo')+1]
    for fi in files:
        if fi.endswith('.npz'):
            # Load file :

            loadname = path.join(st.path, 'database/'+rep, fi)
            print('-> Bipolarize: '+loadname)
            mat = np.load(loadname, allow_pickle=True)
            print('all files', mat.files)

            #~ print(mat['channel'])
            # Bipolarize your data :
            x_bip, chan_bip, label_bip, xyz_bip = bipolarization(mat['x'], [mat['channel'][i][0] for i in range(len(mat['channel']))],
                                                                [mat['label'][i][0] for i in range(len(mat['label']))], xyz=mat['xyz'],)
            #~ print(mat['channel'])
            # Save bipolarized data :
            mat = dict(mat)
            mat['x'], mat['channel'], mat['xyz'], mat['label'] = x_bip, chan_bip, xyz_bip, label_bip
            savename = loadname.replace('.npz', '_bipo.npz')
            np.savez(savename, **mat)
