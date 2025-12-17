"""This script illustrate how you can simply plot your data."""
import numpy as np
import matplotlib.pyplot as plt
from os import path

from brainpipe.system import study
from brainpipe.visual import BorderPlot, addLines
from brainpipe.feature import power

# Define which condition / electrode / frequency to visualize
st = study('Olfacto')
file = 'CHAF_R123_OrdreMFconcat_trigg01_bipo.npz'
elec = 45  # The electrode to plot
# Power settings :
f, fname = [4, 8], 'theta'
baseline = (10, 1536)
norm = 3
# Load file :
loadname = path.join(st.path, 'database', file)
mat = np.load(loadname)
x, sf, channel = mat['x'][elec, ...], int(mat['sf']), mat['channel']
# Get power :
#print('SHAPE OF x: ', x.shape)
powObj = power(sf, x.shape[0], f=f, baseline=baseline, norm=norm)
xpow = np.squeeze(powObj.get(x[np.newaxis, ...])[0])
#xpow = np.squeeze(powObj.get(x)[0])
# Finally plot you power :
time = 1000 * np.arange(x.shape[0]) / sf
BorderPlot(time, xpow, color='slateblue', kind='sem', xlabel='Time (ms)',
           ylabel='Power', title='Power in ' + fname + ' band for electrode '+channel[elec])
addLines(plt.gca(), vLines=baseline, vColor=[
         'firebrick'] * 2, vWidth=[2] * 2)
plt.show()
