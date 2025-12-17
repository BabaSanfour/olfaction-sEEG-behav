# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import timefreqtools

sr = 1000.
times = np.arange(0, 10, 1./sr)
sig = np.sin(2*np.pi*5.*times) * 1./5. + np.sin(2*np.pi*10.*times) * 1./10 + np.sin(2*np.pi*80*times) * 1/80.


f_start, f_stop = 0, 100.
wt, times, freqs, sr2 = timefreqtools.compute_timefreq(sig, sr, f_start, f_stop, delta_freq=0.1,
            f0=10,  normalisation =1, nb_freq=None, min_sampling_rate=1000, returns='all')

spectrum = np.mean(np.abs(wt), axis=0)
fig, ax = plt.subplots()
ax.plot(freqs, spectrum)

fig, ax = plt.subplots()
timefreqtools.plot_tfr(wt, times, freqs, ax=ax, colorbar = True)

plt.show()