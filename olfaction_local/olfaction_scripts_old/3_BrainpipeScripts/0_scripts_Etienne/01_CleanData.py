"""This script will clean your data/channels/coordinates,
adapt it from matlab format and save in a convenient
numpy format.
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


if __name__ == '__main__':
    # Filesettings :
    st = study('Olfacto') # Load the study
    pathdata = r'C:\Users\Anne-Lise\Dropbox\Intra_EM\2_eeg2mat\data_mat\Concat_files' # path where data are located
    trigger = ["01"] # list of triggers
    suj = [('PIRJ', 'MF')] # list of tuples containing (subjects, encoding)

    for tr in trigger:
        for su in suj:
            dico, shortname = _cleanData(pathdata, su[0], su[1], tr)
            path2save = path.join(st.path, 'database', shortname+'.npz') # make sure it's windows compatible
            np.savez(path2save, **dico)
            print('-> Save as : '+path2save+'\n')
