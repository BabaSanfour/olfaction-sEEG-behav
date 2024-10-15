
from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

##############################################################################
roi, freq, exp, pval = 'IFG', 'theta', 'Enc', '0.05'
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
#path_npz = join(path, 'figure/TPSim_LDA_'+exp+'_by_cond_6freqs_3s_dissim/npy_figs/')
path_npz = join(path, 'feature/TPSim_3groups_'+exp+'/Ttest_odor_identity_v=1_elecs=all/')
path_to_save = join(path_npz)
#path_to_save = join(path_npz, 'da_plots/')
#f_form = join(path_npz, '{}_{}_{}_{}_tps.csv'.format(freq,roi,exp,pval))
#f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_da_brain.png')
f_form = join(path_npz, 'Bilan_All_subjects_Ttests_1samp_IFG_bonf_0.05_theta.csv')
#f_form = join(path_npz, 'Bilan_All_subjects_OLS_btw_'+freq+'_mem_groups_'+roi+'_bonf_p0.05.csv')
f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_Tvals.png')
rotation = ['right','bottom']
clim_da = (-0.1,0)
###############################################################################

# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(800, 400))

# =============================================================================
#						Electrodes sources by subject - TOP VIEW
# =============================================================================
# #Define a source object with a different color by subject
df = pd.read_csv(f_form)
x,y,z = df['x'].values[:,np.newaxis], df['y'].values[:,np.newaxis], df['z'].values[:,np.newaxis]
xyz = np.concatenate((x,y,z),axis=1)
#data = df['Tvals'].values
data = df['Tvals_btw'].values
print(data)

# # =============================================================================
# #						Electrodes sources by subject - BOTTOM VIEWS
# # =============================================================================

#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
b_obj_r.alpha = 0.09
s_obj_r = SourceObj('Rep_left', xyz, data=data)
b_obj_r.project_sources(s_obj_r,'modulation', cmap='autumn_r', clim=clim_da,radius=10)
s_obj_r.visible_obj = False
sc.add_to_subplot (b_obj_r, row=0, col=0, rotate=rotation[0])
sc.add_to_subplot(s_obj_r, row=0, col=0)

#Define a brain object to plot
b_obj_2 = BrainObj('B2', translucent=False, hemisphere='both')
b_obj_2.alpha = 0.09
s_obj_2 = SourceObj('modulation', xyz, data=data)
b_obj_2.project_sources(s_obj_2, 'modulation', cmap='autumn_r', clim=clim_da, radius=10)
s_obj_2.visible_obj = False
sc.add_to_subplot(b_obj_2, row=0, col=1, rotate=rotation[1])
sc.add_to_subplot(s_obj_2, row=0, col=1)

sc.screenshot(f_form_save.format('all'),autocrop=True,print_size=(6,7))
