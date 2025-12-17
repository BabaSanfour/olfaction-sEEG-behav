
from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
from utils import rename_elecs
import pandas as pd
from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

step = 'FT'
exp,roi,freq = 'Enc', 'MFG', 'theta'
rois_to_keep = ['ACC','Amg','Amg-PirT','HC','IFG','Ins','MFG','OFC','PHG',
                'SFG','pPirT']
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
#path_npz = join(path, 'figure/TPS_elecs_plots/')
path_npz = join(path, 'feature/TPSim_3groups_'+exp+'/Ttest_odor_identity_v=1_elecs=all/')
#path_npz = join(path, 'database/Encoding_By_Odor/All_elecs_infos_npz/')
path_to_save = path_npz
f_form = join(path_npz, 'Bilan_All_subjects_Ttests_1samp_'+roi+'_fdr_0.05_theta.csv')
#f_form = join(path_npz, 'All_subjects_correl_OFC_olf_unc_0.05_late_theta_Final_mean=False_pl_L.csv')
f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_proj.png')
#f_form_save = join(path_to_save, 'Brain_proj_distribution_{}_elecs'+step+'_OFC.png')
clim = (1,5)
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
#df_sel = df.loc[df['labels'] == 'OFC_olf']
df_sel = df
#id_subj = np.where([su in ['S0','S1','S2','S4','S5','S6'] for su in arch_sig['su_codes'][id_rois]])
#xyz, subjects = arch_sig['s_xyz'][id_rois][id_subj],arch_sig['su_codes'][id_rois][id_subj]
xyz = np.concatenate((df_sel['x'][:,np.newaxis],df_sel['y'][:,np.newaxis],df_sel['z'][:,np.newaxis]),axis=1)
#xyz, subjects = arch_sig['s_xyz'],arch_sig['su_codes']
# b_obj_top = BrainObj('B2', translucent=False, hemisphere='both')
# b_obj_top.alpha = 0.09
# s_obj_top = SourceObj('Rep_right', xyz)
# s_obj_top.set_visible_sources('all')
# s_obj_top.visible_obj = False
# b_obj_top.project_sources(s_obj_top,'repartition', cmap='inferno', clim=(1,10))
# sc.add_to_subplot(b_obj_top, row=0, col=0, rotate='bottom',title='Electrodes repartition all subject')
# sc.add_to_subplot(s_obj_top, row=0, col=0)

# # # =============================================================================
# # #						Electrodes sources by subject - RIGHT VIEW
# # # =============================================================================
# #Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
b_obj_r.alpha = 0.09
s_obj_r = SourceObj('Rep_right', xyz)
#s_obj_r.set_visible_sources('right')
b_obj_r.project_sources(s_obj_r,'repartition', cmap='inferno', clim=clim,radius=10)
s_obj_r.visible_obj = False
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
sc.add_to_subplot(s_obj_r, row=0, col=1)

# # # =============================================================================
# # #						Electrodes sources by subject - LEFT VIEW
# # # =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=False, hemisphere='left')
b_obj_l.alpha = 0.09
s_obj_l = SourceObj('Rep_left', xyz)
s_obj_l.set_visible_sources('left')
b_obj_l.project_sources(s_obj_l,'repartition', cmap='inferno', clim=clim,radius=10)
s_obj_l.visible_obj = False
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
sc.add_to_subplot(s_obj_l, row=0, col=2)
#
# b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='right')
# b_obj_l2.alpha = 0.09
# s_obj_l2 = SourceObj('Rep_left', xyz)
# s_obj_l2.set_visible_sources('left')
# b_obj_l2.project_sources(s_obj_l2,'repartition', cmap='inferno', clim=(1,6),radius=10)
# s_obj_l2.visible_obj = False
# sc.add_to_subplot(b_obj_l2, row=0, col=3, rotate='right')
# sc.add_to_subplot(s_obj_l2, row=0, col=3)
#cb_rep = ColorbarObj (b_obj_l, cblabel='Number of electrodes', txtcolor='black', **CBAR_STATE)
#sc.add_to_subplot(cb_rep, row=0, col=3, width_max=200, height_max=600, )

# # =============================================================================
# #						Electrodes sources by subject - MEDIAL VIEWS
# # =============================================================================
#Define a brain object to plot
# b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
# b_obj_r.alpha = 0.09
# s_obj_r = SourceObj('Rep_right', xyz)
# s_obj_r.set_visible_sources('all')
# b_obj_r.project_sources(s_obj_r,'repartition', cmap='inferno', clim=(1,10))
# s_obj_r.visible_obj = False
# sc.add_to_subplot(b_obj_r, row=0, col=0, rotate='left')
# sc.add_to_subplot(s_obj_r, row=0, col=0)

#Define a brain object to plot
# b_obj_r = BrainObj('B2', translucent=False, hemisphere='left')
# b_obj_r.alpha = 0.09
# s_obj_r = SourceObj('Rep_left', xyz)
# s_obj_r.set_visible_sources('left')
# b_obj_r.project_sources(s_obj_r,'repartition', cmap='inferno', clim=(1,10))
# s_obj_r.visible_obj = False
# sc.add_to_subplot(b_obj_r, row=0, col=3, rotate='right')
# sc.add_to_subplot(s_obj_r, row=0, col=3)

#sc.preview()
sc.screenshot(f_form_save.format('all'),autocrop=True,print_size=(6,7))
