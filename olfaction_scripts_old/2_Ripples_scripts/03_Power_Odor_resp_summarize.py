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
PATH_FEAT = join(st.path, 'feature/Encoding_Retrieval/Wilcox_pre_post_avg/')
feature, corr, thr = 'xpow', 'bonf', 0.05
fnames = ['HFA','ripple']
path_save = PATH_FEAT+'Wilcox_odor_responsive_{}_corr={}_thr={}_roi={}_sign={}.csv'

"""
This script summarize electrodes significant in Wilcoxon non-parametric Ttest
"""

for f,freq in enumerate(fnames):
    print('--> processing {} corr={} thr={}'.format(freq,corr,thr))
    df = pd.read_csv(PATH_FEAT+'Wilcox_odor_responsive_{}_corr={}_thr={}.csv'.format(freq,corr,thr))

    sig_rois = np.unique(df['new_labels'])
    vip_rois = ['olf','OFC_OLF','Ins_OLF','orbital']
    for roi in sig_rois:
        df_roi_f = df.loc[df['new_labels'] == roi]
        df_roi_f['sign'] = np.sign(df_roi_f['Tvals'])
        inc = (df_roi_f.loc[df_roi_f.sign == 1.0]).shape[0]
        dec = (df_roi_f.loc[df_roi_f.sign == -1.0]).shape[0]
        df_inc = df_roi_f.loc[df_roi_f.sign == 1.0].groupby(['subjects']).count()

        if (df_inc.shape[0] >= 3) or (df_inc.shape[0] >=2 and roi in vip_rois):
            print(roi,'%s electrodes showed increased TPSim, while %s showed decrease out of %s elecs'
                  % (inc,dec,inc+dec))
            df_plot = df_roi_f.loc[df_roi_f.sign == 1.0]
            df_plot.to_csv(path_save.format(freq,corr,str(thr),roi,'inc'))
