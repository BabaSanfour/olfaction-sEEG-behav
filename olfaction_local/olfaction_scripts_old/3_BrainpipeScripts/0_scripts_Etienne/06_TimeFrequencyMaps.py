"""This script illustrate how to compute and plot time-frequency maps."""
import numpy as np
import matplotlib.pyplot as plt
from os import path

from brainpipe.system import study
from brainpipe.visual import BorderPlot, addLines
from brainpipe.feature import TF
	
if __name__ == '__main__':
    # Define which condition / electrode / frequency to visualize
    st = study('Olfacto')
    file = 'CHAF_R123_OrdreMFconcat_trigg01_bipo.npz'
    elec = 45  # The electrode to plot
	
    # Power settings :
    f = (0.5, 200, 4, 2)  # Frequency vector: (from, to, width, step)
    baseline = (10, 250)
    norm = 3
    width, step = 50, 10
	
    # Load file :
    loadname = path.join(st.path, 'database', file)
    mat = np.load(loadname)
    x, sf, channel = mat['x'][elec, ...], mat['sf'], mat['channel']
	
    # Get TF :
    time = 1000 * np.arange(x.shape[0]) / sf
    tfObj = TF(sf, x.shape[0], f=f, baseline=baseline,
               norm=norm, time=time, width=width, step=step)
    xpow = np.squeeze(tfObj.get(x)[0])
	
    # Finally plot you power :
    tfObj.plot2D(plt.figure(), 100*xpow, xvec=tfObj.xvec, yvec=tfObj.yvec, xlabel='Time (ms)',
                 ylabel='Frequency (hz)', title='Time-frequency map for channel ' + channel[elec],
                 cblabel='Relative power modulations (%)', cmap='viridis', pltargs={'shading':'gouraud'},
                 vmin=-120, vmax=120, resample=(0.5, 0.1))
    addLines(plt.gca(), vLines=baseline, vColor=[
             'firebrick'] * 2, vWidth=[2] * 2)
    plt.show()
