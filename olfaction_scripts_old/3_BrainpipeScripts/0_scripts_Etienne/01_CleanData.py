"""
This script ingests raw Matlab data and converts it to a clean numpy (.npz) format.

Usage:
    Run this script after setting up 'config.py'.
    It uses 'brainpipe.system.study' to locate the database path.
    Input raw data path is defined in 'config.py'.

Input:
    Matlab files ({subj}_R123_Ordre{encoding}concat_trigg{trigg}.mat) containing:
    - 'x': Raw signal matrix
    - 'chanmonop': Channel info and coordinates
    - 'sf': Sampling frequency

Output:
    Saves .npz files to the 'database/' folder of the Study path.
    The .npz file contains:
    - 'x': Signal (n_channels, n_samples, n_trials)
    - 'xyz': Coordinates
    - 'channel': Channel labels
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from os import path

from brainpipe.system import study

def _cleanData(path, suj, encoding, trigger):
    """The function to clean your data. This is convenient
    because you can loop on your subject/trigger and clean
    everything in one script.
    """
    # Import file :
    file = '{suj}_R123_Ordre{encoding}concat_trigg{trigg}.mat'
    loadname = file.format(suj=suj, encoding=encoding, trigg=trigger)
    print('-> Processing '+loadname)
    mat = loadmat(path.join(pathdata, loadname))

    # Make a clean version of channel name / xyz :
    chan = mat['chanmonop'][0]
    channel = [k[0][0] for k in chan[0]]  # Clean channel name
    xyz = np.zeros((len(channel), 3), dtype=float)
    for k in range(len(channel)):
        xyz[k, 0] = chan[1][k][0]
        xyz[k, 1] = chan[2][k][0]
        xyz[k, 2] = chan[3][k][0]

    # Clean matrix of data (remove unecessary channels)
    x = mat['x']
    if x.ndim == 2:  # This test is in case of x being a 2D matrix (make it 3D)
        x = x[..., np.newaxis]

    # Update variable inside the loaded file :
    dico = dict(mat)
    dico['x'] = x[0:len(channel), ...]
    dico['xyz'] = xyz
    dico['channel'] = channel
    dico['sf'] = dico['sf'][0][0]
    del dico['chanmonop']
    del dico['__header__']
    del dico['__globals__']
    del dico['__version__']

    return dico, loadname.split('.mat')[0]


import sys
from pathlib import Path

if __name__ == '__main__':
    # Add parent directory to path to import config
    current_dir = Path(__file__).resolve().parent
    # olfaction_scripts_old is 2 levels up
    root_dir = current_dir.parent.parent
    sys.path.append(str(root_dir))
    
    try:
        import config
    except ImportError:
        # Create a dummy config if not found to avoid crashing
        class config:
            DEFAULT_RAW_PATH = r'./data'
            STUDY_NAME = 'Olfacto'

    # Filesettings :
    st = study(config.STUDY_NAME) # Load the study
    pathdata = config.get_raw_data_path() # Path where raw Mat files are located
    
    print(f"-> Looking for raw data in: {pathdata}")
    
    trigger = ["01"] # list of triggers
    suj = [('PIRJ', 'MF')] # list of tuples containing (subjects, encoding)

    for tr in trigger:
        for su in suj:
            dico, shortname = _cleanData(pathdata, su[0], su[1], tr)
            path2save = path.join(st.path, 'database', shortname+'.npz') # make sure it's windows compatible
            np.savez(path2save, **dico)
            print('-> Save as : '+path2save+'\n')
