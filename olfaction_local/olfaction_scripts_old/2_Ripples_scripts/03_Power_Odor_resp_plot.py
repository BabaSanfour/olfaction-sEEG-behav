"""
PLOT by ROI results from HFA and Ripple ranges through time
note : T values for all electrodes are projected on 3D MNI brain with Visbrain
"""
import numpy as np
from os.path import join, exists
from os import makedirs, listdir
from itertools import product
import pandas as pd

from brainpipe.system import study
import matplotlib.pyplot as plt
from brainpipe.visual import BorderPlot, addLines, rmaxis
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from utils_functions import compute_big_npz

if __name__ == "__main__":

    st = study('Ripples')
    PATH_FEAT = join(st.path, 'feature/Encoding_Retrieval/')
    PATH_FIG = join(st.path, 'feature/Encoding_Retrieval/Wilcox_pre_post_avg/')
    npz_file = join(PATH_FEAT, 'All_subj_all_odors_feat_Zscore.npz')
    plt_save = join(PATH_FIG, 'Pow_all_subj_feat={}_roi={}_th={}.png')

    feature, corr, thr = 'xpow', 'bonf', 0.05
    fnames = ['HFA', 'ripple']
    st_searchs = ['feat.npz', 'feat_Zscore.npz']
    val_concat = ['xpow', 'channel', 'label', 'xyz']
    val_unique = ['f','split','step','width','fname','time']
    to_plot = [50,100]
    step = 'plot_by_roi' #'compute_npz' ,'plot_by_roi'

    if step == 'plot_by_roi': #load npz and select electrodes to plot
        for fname in fnames:
            files = [f for f in listdir(PATH_FIG) if (fname in f) & ('roi' in f) & (f.endswith('.csv'))]
            rois = [f.split('=')[3][:-5] for f in files]

            for fi,roi in zip(files, rois):
                print('processing -->', fname, roi)
                df_res = pd.read_csv(PATH_FIG+fi)
                mat = np.load(npz_file, allow_pickle=True)

                data_plot = np.array([])
                for su, ch in zip(df_res['subjects'],df_res['channels']):
                    id_f = [i for i,freq in enumerate(mat['fname']) if freq == fname]
                    idx_ch = [i for i,chan in enumerate(mat['channel']) \
                            if (chan == ch) & (mat['subjects'][i] == su)]
                    data = mat['xpow'][id_f,idx_ch,to_plot[0]:to_plot[1]]
                    data_plot = np.concatenate((data_plot,data), axis=0) if np.size(data_plot) else data

                xfmt = ScalarFormatter(useMathText=True)
                xfmt.set_powerlimits((0,3))
                fig = plt.figure(1,figsize=(7,7))
                title = fname+' - Power dynamics in '+roi +'('+str(data_plot.shape[0])+' elecs)'
                fig.suptitle(title, fontsize=12)

                plt.subplot(211)
                time = mat['time'][to_plot[0]:to_plot[1]] - 3
                BorderPlot(time, data_plot, kind='sem',color=['blue','c','c'],linewidth=2,
                              ncol=1, xlabel='Time (s)',ylabel = r'Power',
                              legend=['Odor_m'])
                rmaxis(plt.gca(), ['right', 'top'])
                addLines(plt.gca(), vLines=[0], vColor=['darkgray'], vWidth=[2])
                plt.legend(loc=0, handletextpad=0.1, frameon=False)
                plt.gca().yaxis.set_major_locator(MaxNLocator(3,integer=True))

                plt.subplot(212)
                plt.plot(time, data_plot.swapaxes(0,1))
                #plt.ylim([-1,3])
                plt.legend([su+'_'+ch for su,ch in zip(df_res['subjects'],df_res['channels'])])
                plt.savefig(plt_save.format(fname,roi,thr))
                plt.clf()
                plt.close()

    if step == 'compute_npz': #create big file with all data
        for st_search in st_searchs:
            big_npz = compute_big_npz(PATH_FEAT, st_search, val_concat, val_unique)
            np.savez(PATH_FEAT+'All_subj_all_odors_'+st_search,**big_npz)
