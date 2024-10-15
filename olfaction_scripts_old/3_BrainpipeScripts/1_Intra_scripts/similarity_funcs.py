from os.path import join, exists
from os import makedirs

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations, product

from brainpipe.system import study
from utils import subjects, context_su, odor_groups_wgth, odor_groups_3wgth


def compute_tps_wth(xpow_od):
    nchs, npts, ntrials = xpow_od.shape
    R_elecs, p_elecs, trials_list = np.array([]), np.array([]), []
    for elec in range(nchs):
        pow0 = xpow_od[elec,...]
        R_trials, p_trials = np.array([]), np.array([])
        for t0, t1 in combinations(np.arange(ntrials), 2):
            R, p = stats.pearsonr(pow0[:,t0],pow0[:,t1])
            R_trials = np.vstack((R_trials,R)) if np.size(R_trials) else np.array(R)
            p_trials = np.vstack((p_trials,p)) if np.size(p_trials) else np.array(p)
            trials_list.extend(str(t0)+'_'+str(t1))
        R_elecs = np.vstack((R_elecs,R_trials.T)) if np.size(R_elecs) else R_trials.T
        p_elecs = np.vstack((p_elecs,p_trials.T)) if np.size(p_elecs) else p_trials.T
    return R_elecs, p_elecs, trials_list

def compute_tps_btw(x_od, x_od2):
    R_freq, p_freq = np.array([]), np.array([])
    nelecs, npts, ntrials_od = x_od.shape
    _, _, ntrials_od2 = x_od2.shape
    for elec in range(nelecs):
        pow1, pow2 = x_od[elec], x_od2[elec]
        R_trials, p_trials = np.array([]), np.array([])
        for t0, t1 in product(range(ntrials_od),range(ntrials_od2)):
            R, p = stats.pearsonr(pow1[:,t0],pow2[:,t1])
            R_trials = np.vstack((R_trials,R)) if np.size(R_trials) else np.array(R)
            p_trials = np.vstack((p_trials,p)) if np.size(p_trials) else np.array(p)
        R_freq = np.vstack((R_freq,R_trials.T)) if np.size(R_freq) else R_trials.T
        p_freq = np.vstack((p_freq,p_trials.T)) if np.size(p_freq) else p_trials.T
    return R_freq, p_freq
