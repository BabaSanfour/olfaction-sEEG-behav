"""This script illustrate how to use visbrain to plot your power features on a
standard MNI brain.

Need to install visbrain (from my repo) and paste the folder into
.../site-packages/. You will probably have a problem with pyqt. In that
case, you need to install pyqt 4 (something like this: conda install
pyqt='4.17')

"""
import numpy as np
from os import path

from brainpipe.system import study
from visbrain import vbrain

if __name__ == '__main__':
    # Load one power file :
    st = study('Olfacto')
    pattern = '_expect_E1E2_concat_allfilter1_bipo_meanpow_alpha.npz'  # Pattern for loading all your subjects
    files = st.search(pattern, folder='feature/Power_Encoding_Odor_rest_th40_art400_30_250/')
    ## [2, 4], [5, 7], [8, 13], [13, 30], [30, 60], [60, 120]] Frequency vector
    ## fname = ['delta', 'theta', 'alpha', 'beta', 'gamma30-60', 'gamma60-120']
    frq2plt = 1 # The frequency to plot
    window = 22  # Select one window
    # Now, we construct the data to plot and coordinates :
    s_data = np.array([])
    s_xyz = np.array([])
    for fi in files:
        # Detect the subject :
        suj = fi.split('_E1E2')[0]
        # Load only xyz for this subject :
        xyzFile = st.search(suj, '_bipo', folder='database/TS_E_all_cond_by_block_trigs_th40_art400_30_250/')[0]
        print('-> Coordinates file: ' + xyzFile)
        xyz = np.load(path.join(st.path, 'database/TS_E_all_cond_by_block_trigs_th40_art400_30_250/', xyzFile))['xyz']
        # Now load the power features :
        print('-> Feature file: ' + fi)
        mat = np.load(path.join(st.path, 'feature/Power_Encoding_Odor_rest_th40_art400_30_250/', fi))
        x = mat['xpow'][frq2plt, :, window, :].mean(1)
        fname = mat['fname']
        # Add x and xyz to sources variables :
        s_data = np.concatenate((s_data, x)) if s_data.size else x
        s_xyz = np.concatenate((s_xyz, xyz)) if s_xyz.size else xyz
    # Create a visbrain instance :
    ui_savename = suj + '_Fcy-' + fname[frq2plt] + '_Window-' + str(window)
    vb = vbrain(s_data=s_data, s_xyz=s_xyz, ui_savename=ui_savename)
    vb.show()
