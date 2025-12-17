from itertools import product, combinations
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.sparse import eye


def _compute_rsa(trial1, trial2, stat):
    """Computes the Representational Similarity Analysis (RSA) between two time
        frequency maps.

    Parameters
    ----------
    trial1, trial2 : arrays
        Must be of shape(n_freq x n_samples).
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.

    Returns
    -------
    rsa : array
        The Representational Similarity matrix of shape (n_samples x n_samples).
    """
    corr = pearsonr if stat == "pearson" else spearmanr

    assert (
        trial1.shape == trial2.shape
    ), "Error: shape of trial1 and trial2 must be the same"

    n_samples = trial1.shape[-1]

    rsa = np.zeros((n_samples, n_samples))
    for i, sample1 in enumerate(trial1.T):
        for j, sample2 in enumerate(trial2.T):
            rsa[i, j] = corr(sample1, sample2)[0]

    #reverse the order of the rows to have low origin in rsa matrices
    rsa = rsa[::-1]

    # return Fisher-z transform
    return np.arctanh(rsa)


def _btw_rsa_twoconds(cond1, cond2, stat,average=True):
    """Computes btw rsa of the data.

    Parameters
    ----------
    cond1, cond2 : array
        the data of shape (n_freq x n_samples x n_trials).
        WARNING: ASSUMES cond1 AND cond2 HAVE THE nfreq and nsamples.
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.

    Returns
    -------
    rsa : array
        The Representational Similarity matrix of shape 
        (n_samples x n_samples) if average = True
        OR (n_trials x n_samples x n_samples) if average = False
    """
    assert (
        cond1.shape[:-1] == cond2.shape[:-1]
    ), "Error: cond1 and cond2 should be of the same dimensions for freq and samples"
    n_samples = cond1.shape[1]

    n_comb = len([c for c,_ in product(cond1.T,cond2.T)])
    print('nb of combinations of trials computed', n_comb)
    if average == True:
        rsa = np.zeros((n_samples, n_samples))
        for trial1, trial2 in product(cond1.T, cond2.T):
            rsa += _compute_rsa(trial1.T, trial2.T, stat=stat)
        rsa = rsa / n_comb
    else :
        rsa = np.ones((n_comb, n_samples, n_samples))
        index = 0
        for trial1, trial2 in product(cond1.T,cond2.T):
            rsa[index] = _compute_rsa(trial1.T,trial2.T,stat=stat)[np.newaxis]
            index += 1
    return rsa


def within_rsa(data, stat="pearson"):
    """Computes within rsa of the data.

    Parameters
    ----------
    data : array
        The data of shape (n_freq x n_samples x n_trials).
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.

    Returns
    -------
    rsa : array
        The Representational Similarity matrix of shape (n_samples x n_samples).
    """
    n_samples, n_trials = data.shape[1:]
    rsa = np.zeros((n_samples, n_samples))
    counter = 0
    for i, j in combinations(np.arange(n_trials), 2):
        rsa += _compute_rsa(data[:, :, i], data[:, :, j], stat=stat)
        counter += 1
    return rsa / counter


def btw_rsa(conds, stat, rep=1, n_jobs=-1):
    """

    Parameters
    ----------
    conds : list of arrays
        Each array must be of shape (n_freq x n_samples x n_trials).
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.
    rep : int
        The number of bootstrapping steps. Use if data is unbalanced.


    Returns
    -------
    rsa : array
        The Representational Similarity matrices of shape (n_channels x n_samples x n_samples).

    """
    assert len(conds.shape) == 3, "Error: conds should have only 3 dimensions"
    n_freq, n_samples = conds[0].shape[:-1]

    # Do we need bootstrapping ?
    for _ in range(rep):
        # get the minimum number of stim in our conditions
        mini = min([cond.shape[-1] for cond in conds])

        # preallocation
        balanced_conds = [np.zeros((n_freq, n_samples, mini)) for _ in conds]

        # selecting random trials from each condition and concatenating in one matrix
        for i, cond in enumerate(conds):
            index = np.random.choice(np.arange(cond.shape[-1]), mini, replace=False)
            balanced_conds[i] = conds[i][:, :, index]

        rsa = np.zeros((n_samples, n_samples))
        counter = 0
        for i, j in combinations(np.arange(len(conds)), 2):
            rsa += _btw_rsa_twoconds(balanced_conds[i], balanced_conds[j], stat=stat)
            counter += 1
    return rsa / counter


if __name__ == "__main__":
    N_BOOT = 2

    N_CHANNELS = 20
    N_SUBJ = 16
    N_FREQ = 5
    N_SAMPLES = 10
    N_TRIALS = 10
    N_TRIALS2 = 5
    N_TRIALS3 = 3
    # n_freq x n_channels x n_samples x n_trials
    cond1 = np.random.rand(N_FREQ, N_SAMPLES, N_TRIALS)
    cond2 = np.random.rand(N_FREQ, N_SAMPLES, N_TRIALS2)
    cond3 = np.random.rand(N_FREQ, N_SAMPLES, N_TRIALS2)
    conds = [cond1, cond2, cond3]
    for stat in ("spearman", "pearson"):
        print("testing", stat)
        btw_rsa(conds, stat)[np.newaxis]
        within_rsa(cond1, stat)[np.newaxis]

    from mne.stats import permutation_cluster_test

    rsa1 = np.array(
        [
            [
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 8, 6, 0, 0, 0, 0, 0],
                [0, 4, 2, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
            ]
            for _ in range(N_SUBJ)
        ]
    )
    rsa2 = np.array(
        [
            [
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 1, 0, 0, 0, 0, 0],
                [0, 1, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
            ]
            for _ in range(N_SUBJ)
        ]
    )

    # connectivity = eye(8, 8, 0) + eye(8, 8, 1) + eye(8, 8, -1)
    connectivity = sum([eye(8, 8, i) for i in range(-1, 2)])
    # connectivity = None
    _, clu, p, _ = permutation_cluster_test([rsa1, rsa2], connectivity=connectivity)
    print(len(clu), clu[0].shape, clu[0], sep="\n")
