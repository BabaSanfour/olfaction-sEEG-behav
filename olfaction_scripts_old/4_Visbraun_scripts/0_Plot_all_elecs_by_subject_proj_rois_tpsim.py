
from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj, RoiObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
roi, freq, exp, pval = 'aHC', 'theta', 'Enc', '0.05'
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'feature/TPSim_3groups_'+exp+'/Ttest_odor_identity_v=1_elecs=all/')
#path_npz = join(path, 'figure/TPSim_LDA_'+exp+'_by_cond_6freqs_3s_dissim/npy_figs/')
path_to_save = join(path_npz)
#f_form = join(path_npz, 'TPS_self_Enc_'+freq+'_'+roi+'_elecs_sig_coords.npy')
f_form = join(path_npz, 'Bilan_All_subjects_Ttests_1samp_aHC_fdr_0.05_theta.csv')
#data_form = join(path_npz, 'TPS_self_Enc_'+freq+'_'+roi+'_pvals.npy')
f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_Tvals.png')
#cb_save = join(path_to_save, 'Rep_Signif_all_Brain_CB_EpiScore.png')
roi_plot = 'Hippocampus'#'Hippocampus'#'INS_v_PIsul'
clim = (0,0.1)#(0,0.8)
##############################################################################

# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(1100, 400))

# =============================================================================
#						Electrodes sources by subject - TOP VIEW
# =============================================================================
# #Define a source object with a different color by subject
df = pd.read_csv(f_form)
x,y,z = df['x'].values[:,np.newaxis], df['y'].values[:,np.newaxis], df['z'].values[:,np.newaxis]
xyz = np.concatenate((x,y,z),axis=1)
#data = df[freq+'_AUC'].values
#data = df['Tvals'].values
data = df['Tvals_btw'].values
print(data)

b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
b_obj_top.alpha = 0.09
sc.add_to_subplot(b_obj_top, row=0, col=0, rotate='front',title='LTM electrodes repartition')
s_obj_top = SourceObj('Rep LTM', xyz, data=data)
s_obj_top.visible_obj = False
sc.add_to_subplot(s_obj_top, row=0, col=0)

roi_obj = RoiObj('aal')
idx = roi_obj.where_is(roi_plot)
#roi_obj.select_roi(select=dmn_idx)
#roi_obj = RoiObj('aal', border=False)
#idx = roi_obj.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj.select_roi(select=idx, unique_color=False, smooth=2)
roi_obj.project_sources(s_obj_top, 'modulation', cmap='autumn_r', clim=clim)
sc.add_to_subplot(roi_obj, row=0, col=0)
#sc.preview()

# # =============================================================================
# #						Electrodes sources by subject - RIGHT VIEW
# # =============================================================================
#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
b_obj_r.alpha = 0.09
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
s_obj_r = SourceObj('Rep_right', xyz, data=data)
s_obj_r.set_visible_sources('right')
s_obj_r.visible_obj = False
sc.add_to_subplot(s_obj_r, row=0, col=1)
roi_obj_r = RoiObj('aal', border=False)
idx = roi_obj_r.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj_r.select_roi(select=idx, unique_color=False, smooth=2)
roi_obj_r.project_sources(s_obj_r, 'modulation', cmap='autumn_r', clim=clim)
sc.add_to_subplot(roi_obj_r, row=0, col=1)

# # =============================================================================
# #						Electrodes sources by subject - LEFT VIEW
# # =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
b_obj_l.alpha = 0.09
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
s_obj_l = SourceObj('Rep_left', xyz, data=data)
s_obj_l.set_visible_sources('left')
s_obj_l.visible_obj = False
sc.add_to_subplot(s_obj_l, row=0, col=2)
roi_obj_l = RoiObj('aal', border=False)
idx = roi_obj_l.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj_l.select_roi(select=idx, unique_color=False, smooth=2)
roi_obj_l.project_sources(s_obj_l, 'modulation', cmap='autumn_r', clim=clim)
sc.add_to_subplot(roi_obj_l, row=0, col=2)

cb_rep = ColorbarObj(roi_obj_l, cblabel='Number of electrodes', **CBAR_STATE)
sc.add_to_subplot(cb_rep, row=0, col=3, width_max=200, height_max=600)

#sc.preview()
sc.screenshot(f_form_save,autocrop=True,print_size=(6,7))
