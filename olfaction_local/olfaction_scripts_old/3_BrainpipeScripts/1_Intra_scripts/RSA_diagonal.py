from itertools import product, combinations
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.sparse import eye


def _compute_rsa_diag(trial1, trial2, stat):
    """Computes the Diagonal Similarity Analysis (RSA) between two time
        frequency maps from the same condition

    Parameters
    ----------
    trial1, trial2 : arrays
        Must be of shape(n_freq x n_samples).
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.

    Returns
    -------
    rsa : array
        The Representational Similarity matrix of shape (n_samples).
    """
    corr = pearsonr if stat == "pearson" else spearmanr

    assert (
        trial1.shape == trial2.shape
    ), "Error: shape of trial1 and trial2 must be the same"

    n_samples = trial1.shape[-1]

    rsa = np.zeros((n_samples))
    for i, sample1 in enumerate(trial1.T):
        for sample2 in trial2.T:
            rsa[i] = corr(sample1, sample2)[0]

    # return Fisher-z transform
    return np.arctanh(rsa)


def btw_rsa_twoconds(cond1, cond2, stat):
    """Computes btw rsa of the data.

    Parameters
    ----------
    cond1, cond2 : array
        the data of shape (n_freq x n_samples x n_trials).
        WARNING: ASSUMES cond1 AND cond2 HAVE THE SAME nfreq and nsamples.
    stat : string
        The stat of correlation method used. Can be 'pearson' or 'spearman'.

    Returns
    -------
    rsa : array
        The Representational Similarity matrix of shape 
        (n_combinations x n_samples)
    """
    assert (
        cond1.shape[:-1] == cond2.shape[:-1]
    ), "Error: cond1 and cond2 should be of the same dimensions for freq and samples"
    n_samples = cond1.shape[1]

    n_comb = len([c for c,_ in product(cond1.T,cond2.T)])
    print('nb of combinations of trials computed', n_comb)
    rsa = np.ones((n_comb, n_samples))
    index = 0
    for trial1, trial2 in product(cond1.T,cond2.T):
        rsa[index] = _compute_rsa_diag(trial1.T,trial2.T,stat=stat)
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
        The Representational Similarity matrix of shape (ncombinations, n_samples).
    """
    n_samples, n_trials = data.shape[1:]
    n_comb = len([c for c,_ in combinations(np.arange(n_trials), 2)])

    rsa = np.zeros((n_comb,n_samples))
    index = 0
    for i, j in combinations(np.arange(n_trials), 2):
        rsa[index] += _compute_rsa_diag(data[:, :, i], data[:, :, j], stat=stat)
        index += 1
    return rsa


# def btw_rsa(conds, stat, boot=True, rep=10, n_jobs=-1):
#     """

#     Parameters
#     ----------
#     conds : list of the two arrays
#         Each array must be of shape (n_freq x n_samples x n_trials).
#     stat : string
#         The stat of correlation method used. Can be 'pearson' or 'spearman'.
#     boot : bool
#         Bootstrapping rsa computations when conditions are unbalanced.
#     rep : int
#         The number of bootstrapping steps when bootstrap is True.


#     Returns
#     -------
#     rsa : array
#         The Representational Similarity matrices of shape (n_samples).

#     """
#     assert len(conds[0].shape) == 3, "Error: conds should have only 3 dimensions"
#     n_freq, n_samples = conds[0].shape[:-1]
#     n_trials = [cond.shape[-1] for cond in conds]

#     # Do we need bootstrapping ?
#     if boot == True and n_trials[0]!=n_trials[1]:

#         # get the minimum number of stim in our conditions
#         mini = min([cond.shape[-1] for cond in conds])

#         # preallocation
#         rsa = np.zeros(len(product(np.arange(mini),np.arange(mini))))

#         for counter in range(rep):
#             # preallocation
#             balanced_conds = [np.zeros((n_freq, n_samples, mini)) for _ in conds]

#             # selecting random trials from each condition and concatenating in one matrix
#             for i, cond in enumerate(conds):
#                 index = np.random.choice(np.arange(cond.shape[-1]), mini, replace=False)
#                 balanced_conds[i] = conds[i][:, :, index]

#             rsa += _btw_rsa_twoconds(balanced_conds[0], balanced_conds[1], stat=stat)
#         rsa = rsa / rep

#     else:
#         # preallocation
#         rsa = np.zeros(len(product(np.arange(n_trials[0]),np.arange(n_trials[1]))))
#         rsa += _btw_rsa_twoconds(conds[0], conds[1], stat=stat)

    
#     return rsa


if __name__ == "__main__":
    N_BOOT = 2

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
        #within_rsa(cond1, stat)[np.newaxis]
