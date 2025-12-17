
from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj, RoiObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
exp,roi,freq = 'Enc', 'Amg_pPirT', 'theta'
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'feature/TPSim_3groups_'+exp+'/Ttest_odor_identity_v=1_elecs=all/')
#path_npz = join(path, 'database/Encoding_By_Odor/All_elecs_infos_npz/')
path_to_save = path_npz
#f_form = join(path_npz, '0_all_subjects_info_elecs_sel.csv')
f_form = join(path_npz, 'Bilan_All_subjects_Ttests_1samp_'+roi+'_fdr_0.05_theta.csv')
f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_no_brain_proj.png')
#f_form_save = join(path_to_save, 'Brain_proj_distribution_LTM.png')
###############################################################################
roi_plot, clim, smooth, radius = 'Amygdala', (0,5), 5, 15
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
df_sel = pd.read_csv(f_form)
#df_sel = df.loc[df['labels'] == roi]
xyz = np.concatenate((df_sel['x'][:,np.newaxis],df_sel['y'][:,np.newaxis],df_sel['z'][:,np.newaxis]),axis=1)

b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
b_obj_top.alpha = 0.
sc.add_to_subplot(b_obj_top, row=0, col=0, rotate='front',title='LTM electrodes repartition')
s_obj_top = SourceObj('Rep LTM', xyz)
s_obj_top.visible_obj = False
sc.add_to_subplot(s_obj_top, row=0, col=0)
roi_obj = RoiObj('aal', border=False)
idx = roi_obj.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj.select_roi(select=idx, unique_color=False, smooth=smooth)
roi_obj.project_sources(s_obj_top, 'repartition', cmap='inferno', clim=clim,radius=radius)
sc.add_to_subplot(roi_obj, row=0, col=0)

# # =============================================================================
# #						Electrodes sources by subject - RIGHT VIEW
# # =============================================================================
#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
b_obj_r.alpha = 0.
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
s_obj_r = SourceObj('Rep_right', xyz)
s_obj_r.set_visible_sources('right')
s_obj_r.visible_obj = False
sc.add_to_subplot(s_obj_r, row=0, col=1)
roi_obj_r = RoiObj('aal', border=False)
idx = roi_obj_r.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj_r.select_roi(select=idx, unique_color=False, smooth=smooth)
roi_obj_r.project_sources(s_obj_r, 'repartition', cmap='inferno', clim=clim,radius=radius)
sc.add_to_subplot(roi_obj_r, row=0, col=1)

# # =============================================================================
# #						Electrodes sources by subject - LEFT VIEW
# # =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
b_obj_l.alpha = 0
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
s_obj_l = SourceObj('Rep_left', xyz)
s_obj_l.set_visible_sources('left')
s_obj_l.visible_obj = False
sc.add_to_subplot(s_obj_l, row=0, col=2)
roi_obj_l = RoiObj('aal', border=False)
idx = roi_obj_l.where_is([roi_plot]) #'Amygdala','ParaHippocampal'
roi_obj_l.select_roi(select=idx, unique_color=False, smooth=smooth)
roi_obj_l.project_sources(s_obj_l, 'repartition', cmap='inferno', clim=clim,radius=radius)
sc.add_to_subplot(roi_obj_l, row=0, col=2)

cb_rep = ColorbarObj(roi_obj_l, cblabel='Number of electrodes', **CBAR_STATE)
sc.add_to_subplot(cb_rep, row=0, col=3, width_max=200, height_max=600)

#sc.preview()
sc.screenshot(f_form_save.format('FT'),autocrop=True,print_size=(6,7))
