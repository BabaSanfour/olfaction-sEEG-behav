import numpy as np
from brainpipe.feature import sigfilt
from brainpipe.tools import binarize


def rectification(x, sf, f=[5, 7], width=None, step=None):
    """Computes rectification followed by moving average.
    1. Narrow-band filtering
    2. Absolute value
    3. Moving average    Parameters
    ----------
    x : array of shape (n_electrodes, n_times, n_trials)
    sf : sampling frequency
    f : frequency bands
    width / step : width and step (in samples) of the sliding window

    Returns
    -------
    x_filt : (n_freqs, n_electrodes, n_times, n_trials)
    """
    n_elec, n_times, n_trials = x.shape
    # filter the signal
    s_obj = sigfilt(sf, n_times, f=f)
    x_filt, _ = s_obj.get(x)
    n_freqs = x_filt.shape[0]
    # takes the absolute value
    x_filt = np.abs(x_filt)
    # moving average + down-sample
    if isinstance(width, int) and isinstance(step, int):
        win = binarize(0, n_times, width, step)
        n_win = len(win)
        x_ma = np.zeros((n_freqs, n_elec, n_win, n_trials))
        for n_w, (s, e) in enumerate(win):
            x_ma[:, :, n_w, :] = x_filt[:, :, s:e, :].mean(2)
        return x_ma
    else:
        return x_filt
