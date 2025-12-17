"""This script will clean your data/channels/coordinates,
adapt it from matlab format and save in a convenient
numpy format.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from os import path, listdir, makedirs
from brainpipe.system import study
from params import subjects

def cleanData(rep, filename):
    """The function to clean your data. This is convenient
    because you can loop on your subject/trigger and clean
    everything in one script.
    """

    # Import file :
    print('-> Processing '+filename)
    mat = loadmat(path.join(rep+filename))

    # Make a clean version of channel name / xyz :
    chan_xy = mat['chanmonop'][1:4]
    print (chan_xy.shape)
    chan_names = mat['chanmonop'][0]
    chan_label = mat['chanmonop'][4]
    # channel = [k[0][0] for k in chan[0]]  # Clean channel name
    xyz = np.zeros([len(chan_names), 3], dtype=float)
    for k in range(len(chan_names)):
        xyz[k, 0] = float(chan_xy[0,k][0])
        xyz[k, 1] = float(chan_xy[1,k][0])
        xyz[k, 2] = float(chan_xy[2,k][0])
    xyz = xyz.round(decimals=1, out=None)

    # Clean matrix of data (remove unecessary channels)
    x = mat['x']
    print ('shape of the data', x.shape)
    if x.ndim == 2:  # This test is in case of x being a 2D matrix (make it 3D)
        x = x[..., np.newaxis]
    x = x.swapaxes(0, 1).swapaxes(1,2)
    print (x.shape)

    # Update variable inside the loaded file :
    dico = dict(mat)
    dico['x'] = x[0:len(chan_names), ...]
    dico['xyz'] = (xyz)
    dico['channel'] = chan_names
    dico['label'] = chan_label
    dico['sf'] = dico['sf'][0][0]
    del dico['chanmonop']
    del dico['__header__']
    del dico['__globals__']
    del dico['__version__']

    return dico


if __name__ == '__main__':
    st = study('Ripples')
    step = 'import_data'

    PATH = r'/media/1_Analyses_Intra_EM_Odor/OE_Matrices_NoArt/'
    reps = ['TS_E_all_cond_by_block_trigs_odors/',
            'TS_R_all_cond_by_block_trigs_odors/',
            'TS_R_all_cond_by_block_trigs_rec/']#, 'TS_R_all_cond_by_block_trigs_rec_time/']
    subjects = ['PIRJ','CHAF','VACJ','SEMC','LEFC','FERJ']

    if step == 'import_data':
        for rep in reps:
            files = listdir(PATH+rep+'_mat_by_odor/')
            for f in files:
                dico = cleanData(PATH+rep+'_mat_by_odor/', f)
                path_save = path.join(st.path, 'database/'+rep.split('_')[1]+'_'+rep.split('_')[-1]+'/')
                if not path.exists(path_save):
                    makedirs(path_save)
                np.savez(path_save+f.replace('.mat','.npz'),**dico)
