# -*- coding: utf-8 -*-


import numpy as np
from matplotlib import pyplot
import neo
from OpenElectrophy.timefrequency import TimeFreq
import quantities as pq
import scipy.signal


# construct signal
sr = 500.
times = np.arange(0,4,1./sr, dtype = 'float64')
l = times.size
f1 =1
f2 = 25
f3 = 70
sig1 = np.sin(times*2*np.pi*f1) * 1 * np.hamming(times.size)
sig2 = np.sin(times*2*np.pi*f2)*.3 * np.concatenate([np.zeros(l/2), np.hamming(times.size/4),np.zeros(l/4)])
sig3 = np.sin(times*2*np.pi*f3)*.2 * np.concatenate([np.zeros(l/4), np.hamming(times.size/2),np.zeros(l/4)])
sig = sig1+sig2+sig3 + np.random.randn(l)*.1

neosig = neo.AnalogSignal(sig, units = 'V', sampling_rate = sr*pq.Hz)

# Morlet scalogram
tfr = TimeFreq(neosig, f0 =4.5, f_start=0., f_stop=90, sampling_rate=500.)

# Hilbert
bands = [  (0.1,10.),  (15.,35.), (35,100)]
filtered_sigs = [ ]
hilbert_envelops = [ ]
for f1, f2 in bands:
    Wn = [f1/(sr/2.), f2/(sr/2.) ]
    b, a = scipy.signal.iirfilter(N=3, Wn=Wn, btype = 'bandpass', analog = 0, ftype = 'butter', output = 'ba')
    sigf = scipy.signal.filtfilt(b, a, sig)
    filtered_sigs.append(sigf)
    env = np.abs(scipy.signal.hilbert(sigf))# estimation with hilbert
    hilbert_envelops.append(env)

# Enveloppe with TFR
tfr_envelops = [ ]
for f1, f2 in bands:
    maskf = (tfr.freqs.magnitude>=f1) & (tfr.freqs.magnitude<=f2)
    env = np.max(abs(tfr.map[:,maskf]), axis = 1)
    tfr_envelops.append(env)


#Ploting
fig, axs = pyplot.subplots(nrows=5, sharex = True)

axs[0].plot(times, sig)
tfr.plot(ax = axs[1], colorbar = False)

for sigf in filtered_sigs:
    axs[2].plot(times, sigf)

for env in hilbert_envelops:
    axs[3].plot(times, env)

for env in tfr_envelops:
    print env.shape
    print tfr.times.magnitude.shape
    axs[4].plot(tfr.times.magnitude, env)


pyplot.show()





