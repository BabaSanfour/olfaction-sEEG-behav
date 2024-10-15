"""This script illustrate how to compute and plot time-frequency maps."""
import numpy as np
import matplotlib.pyplot as plt
from os import path
#%matplotlib notebook

from brainpipe.system import study
from brainpipe.visual import BorderPlot, addLines
from brainpipe.feature import TF


# Define which condition / electrode
st = study('Olfacto')
files = st.search('_bipo.npz', folder='database')
#file = 'PIRJ_R123_OrdreMFconcat_trigg00_bipo.npz' # CHANGE fname at the end of the script!!
n_elec = 111
#elec = 10

# Power settings :
f = (0.5, 200, 4, 2)  # Frequency vector: (from, to, width, step)
baseline = (10, 3000)
norm = 3
width, step = 50, 10
window = (3500, 10000) # IN MILISECONDS
trigger = [4000]
for file in files:
    for elec in range(0, n_elec, 1):
        # Load file :
        loadname = path.join(st.path, 'database', file)
        mat = np.load(loadname)
        x, sf, channel = mat['x'][elec, ...], float(mat['sf']), mat['channel']

        # Get TF :
        time = 1000 * np.arange(x.shape[0]) / sf
        tfObj = TF(sf, x.shape[0], f=f, baseline=baseline,
                   norm=norm, time=time, width=width, step=step)
        xtf = np.squeeze(tfObj.get(x)[0])

        # Plot everything that is inside the window parameter :
        timebin = np.array(tfObj.xvec)
        sl = slice(np.argmin(np.abs(timebin-window[0])), np.argmin(np.abs(timebin-window[1])))

        # Finally plot you power :
        tfObj.plot2D(plt.figure(), 100*xtf[:, sl], xvec=tfObj.xvec[sl], yvec=tfObj.yvec, xlabel='Time (ms)',
                     ylabel='Frequency (hz)', title='Time-frequency map for channel ' + channel[elec],
                     cblabel='Relative power modulations (%)', cmap='viridis', pltargs={'shading':'gouraud'},
                     vmin=-120, vmax=120, resample=(0.5, 0.1))
        addLines(plt.gca(), vLines=trigger, vColor=[
                 'firebrick'], vWidth=[2])

        #Save all your plots
        rep = path.join(st.path, 'TimeFrequency')
        #fname = (rep + '\PIRJ_R123_OrdreMFconcat_trigg00_bipo_TF_' + channel [elec] +'.png')
        fname = loadname.replace('.npz', '_TF_'+ channel [elec] +'.png').replace('database', 'TimeFrequency')
        print (fname)
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
        del x, sf, channel
        