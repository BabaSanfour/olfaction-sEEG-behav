import numpy as np
from os.path import join, exists
from os import makedirs

from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt
from brainpipe.feat.utils._feat import _manageWindow
from params import subjects
from utils_functions import rename_elecs_df

from scipy.stats import wilcoxon
from mne.stats import bonferroni_correction, fdr_correction
import pandas as pd

st = study('Ripples')
PATH_FEAT = join(st.path, 'feature/Encoding/')
PATH_SAVE = join(PATH_FEAT, 'Wilcox_pre_post_avg/')
if not exists(PATH_SAVE):
    makedirs(PATH_SAVE)
feature, corr, thr = 'xpow', 'fdr', 0.05
fnames = ['theta']

"""
This script compare pre and post average power using Wilcoxon non-parametric Ttest
Results are computed by freq band and corrected across all electrodes (all subj together)
"""

def create_stats_df(T,p,corr,thr,chans,labs,subjs,freq,x,y,z,all_pre, all_post):

    #create df ALL elecs
    chans_sig, su_sig, lab_sig = np.array(chans), np.array(subjs), np.array(labs)
    Tvals, pvals, x, y, z = np.array(T), np.array(p), np.array(x), np.array(y), np.array(z)
    all_pre, all_post = np.array(all_pre), np.array(all_post)
    df = pd.DataFrame({'subjects':su_sig, 'channels':chans_sig, 'x':x, 'y':y,'z':z,
                'labels':lab_sig,'Tvals':Tvals, 'pvals':pvals, 'pre':all_pre, 'post':all_post})
    df_new = rename_elecs_df(df)
    _, pvals_corr = bonferroni_correction(df_new['pvals'].values,thr) if corr == 'bonf' \
                                                    else fdr_correction(df_new['pvals'].values,thr)
    df_new['p_corr'] = pvals_corr
    df_not_corr = df_new.loc[df_new['pvals']<0.05]
    print('total nb of sig elec', freq, df_not_corr.shape[0], 'on ', df_new.shape[0],'in total')
    df_new.to_csv(PATH_SAVE+'Wilcox_odor_responsive_ALL_{}_corr={}_thr={}.csv'.format(freq,
                        corr,str(thr)),index=False)

    #create signif df
    df_pos = df_new.loc[(df_new['p_corr']<thr)&(df_new['Tvals']>0)]
    df_neg = df_new.loc[(df_new['p_corr']<thr)&(df_new['Tvals']<0)]
    print('total nb of sig elec', freq, 'pos',df_pos.shape[0], 'neg', df_neg.shape[0])
    print(df_pos.groupby(['new_labels','subjects']).sum())
    0/0
    df_sig.to_csv(PATH_SAVE+'/Wilcox_odor_responsive_{}_corr={}_thr={}.csv'.format(freq,
                            corr,str(thr)),index=False)

if __name__ == "__main__":

    for freq in fnames:
        all_Tvals, all_pvals = [], []
        all_chans, all_labels, subj_list = [], [], []
        x, y, z, all_pre, all_post = [], [], [], [], []
        for su in subjects:
            mat = np.load(PATH_FEAT+su+'_cond=ALL_odors=ALL_bipo_INFO_pow.npz', allow_pickle=True)
            id_f = [i for i,f in enumerate(mat['fname']) if f==freq][0]
            data = mat[feature][id_f,...]
            nchans, npts, ntrials = data.shape
            su_list = [su]*nchans
            channels, labels, xyz = mat['channels'], mat['labels'], mat['xyz']
            # baseline -1000 to -300 ms // odor 0 to +2000ms (check mat['time'])
            #print(mat['time'][17:24], mat['time'][27:47])
            pre, post = np.mean(data[:,17:24,:], axis=1), np.mean(data[:,27:47,:],axis=1)
            #print(np.mean(pre,axis=1).shape, post.shape, channels.shape)
            all_pre.extend(np.mean(pre,axis=1)), all_post.extend(np.mean(post,axis=1))

            T_su, p_su = [], []
            for n in range(nchans):
                T, p = wilcoxon(pre[n],post[n])
                T_su.append(T), p_su.append(p)
            all_Tvals.extend(T_su), all_pvals.extend(p_su)
            all_chans.extend(channels), all_labels.extend(labels), subj_list.extend(su_list)
            x.extend(xyz[:,0]), y.extend(xyz[:,1]), z.extend(xyz[:,2])

        print('processing {} using {} correction'.format(freq, corr))
        create_stats_df(all_Tvals, all_pvals, corr,thr,all_chans,
                        all_labels,subj_list,freq,x,y,z, all_pre, all_post)
