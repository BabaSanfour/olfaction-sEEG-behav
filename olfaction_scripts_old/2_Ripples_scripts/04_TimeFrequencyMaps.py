"""This script illustrates how to compute and plot time-frequency maps."""
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, exists
from os import makedirs
import pandas as pd

from brainpipe.system import study
from brainpipe.visual import BorderPlot, addLines
from brainpipe.feature import TF
from utils import subjects

st = study('Ripples')
PATH_FEAT = join(st.path, 'feature/Encoding_Retrieval/')
PATH_DATA = join(st.path, 'database/Encoding_Retrieval/')

def compute_TFmap(su,data,chans,labs):
    """
    Compute TimeFrequency maps features (nfreqs x ntimes)
    Compute statistics with baseline (data normalzed or not)
    Plot and save maps

    Parameters
    ----------
    su : string (name of the subject)
    data : array (nelecs, npts, trials)
    chans,labs : arrays (nelecs)
        channels names and labels of electrodes
    """

    ###############################################################################
    if not exists(PATH_FEAT+'TF_plots/'):
        makedirs(PATH_FEAT+'TF_plots/')
    ###############################################################################

    # Power settings :
    f = (0.5, 120, 4, 2)  # Frequency vector: (from, to, width, step)
    baseline = (1500, 2500)
    window_to_plot = [-1000, 2000] # In seconds
    norm, sf = 3, 500
    width, step = 50, 10

    # Load file :
    for elec in range(data.shape[0]):
        x, chan, lab = data[elec, ...], chans[elec], labs[elec]

        # Get TF :
        time = (1000 * np.arange(x.shape[0]) / sf)-3000
        tfObj = TF(sf, x.shape[0], f=f, baseline=baseline,
                   norm=norm, time=time, width=width, step=step)
        xpow = np.squeeze(tfObj.get(x,n_jobs=-1)[0])

        # Finally plot you power :
        timebin = np.array(tfObj.xvec)
        sl = slice(np.argmin(np.abs(timebin-window_to_plot[0])), np.argmin(np.abs(timebin-window_to_plot[1])))
        title = '{} elec({}) label={} Encoding and Retrieval'.format(su,chan,lab)
        tfObj.plot2D(plt.figure(), 100*xpow[:,sl], xvec=tfObj.xvec[sl], yvec=tfObj.yvec, xlabel='Time (ms)',
                     ylabel='Frequency (hz)', title=title,
                     cblabel='Relative power modulations (%)', cmap='viridis', pltargs={'shading':'gouraud'},
                     vmin=-50, vmax=50, resample=(0.5, 0.1))
        addLines(plt.gca(), vLines=[0], vColor=['firebrick'], vWidth=[2])
        plt.savefig(PATH_FEAT+'TF_plots/'+'{}_TFmap_elec({})_{}.png'.format(su,chan,lab), dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()


if __name__ == "__main__":

    #params to select electrodes to compute
    corr, thr, freq = 'fdr', '0.05', 'gamma_high'
    for su in subjects:
        mat = np.load(PATH_DATA+su+'_all_odors_bip.npz',allow_pickle=True)
        x, labels, chans = mat['x'], mat['label'], mat['channel']
        df = pd.read_csv(PATH_FEAT+'Wilcox_odor_responsive_{}_corr={}_thr={}.csv'.format(freq,corr,thr))
        chan_sig = df['channels'].loc[df['subject']==su]
        id_sel = [i for i,chan in enumerate(chans) if chan in chan_sig.values]
        data_sel, chans_sel, labels_sel = x[id_sel,...], chans[id_sel], labels[id_sel]
        compute_TFmap(su,data_sel,chans_sel,labels_sel)
