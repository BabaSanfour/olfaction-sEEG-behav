from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

step = 'FT'
rois_to_keep = ['ACC','Amg','Amg-PirT','HC','IFG','Ins','MFG','OFC','PHG',
                'SFG','pPirT']
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/figure/'
path_npz = join(path, 'TPS_elecs_plots/')
path_to_save = path_npz
f_form = join(path_npz, 'All_subjects_sources_TPS_elecs.npz')
f_form_save = join(path_npz, 'Brain_distribution_sources_elecs.png')
###############################################################################

# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1,1,1),size=(1100, 400))

# =============================================================================
#						Electrodes sources by subject - TOP VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_top = BrainObj('B2', translucent=True)
b_obj_top.alpha = 0.09
sc.add_to_subplot(b_obj_top, row=0, col=0, 
					title='Electrodes localization by subject - Top')
#Define a source object with a different color by subject
arch_sig = np.load(f_form)
#id_rois = np.where([roi in rois_to_keep for roi in arch_sig['labels']])
#id_subj = np.where([su in ['S0','S1','S2','S4','S5','S6'] for su in arch_sig['su_codes'][id_rois]])
xyz, subjects = arch_sig['s_xyz'],arch_sig['su_codes']
print('number of electrodes:',xyz.shape)
u_color = ["darkblue", "royalblue", "deepskyblue", "yellow", "darkorange","red"]
color = [u_color[int(k[1])] for k in sorted(subjects)]
data = np.arange(len(subjects))
labels = None
#labels = arch_sig['su_codes']#text=labels,text_size=12.,
s_obj_c = SourceObj ('Color', xyz, symbol='disc', color=color,radius_min=10., alpha=.95, data=data,
					text=labels, text_size=12, text_color='white' )
sc.add_to_subplot(s_obj_c, row=0, col=0)
#sc.preview()

# =============================================================================
#						Electrodes sources by subject - RIGHT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
b_obj_r.alpha = 0.09
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right',
					title='Electrodes localization by subject - Right')
s_obj_r = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
			alpha=.95, data=data, text=labels)
s_obj_r.set_visible_sources('right')
# df2 = s_obj_r.analyse_sources(roi_obj='aal',distance=15.,keep_only=['Amygdala (R)', 'Cingulum Ant (R)', 
# 	'Cingulum Mid (R)', 'Frontal Inf Orb (R)', 'Frontal Inf Tri (R)', 'Frontal Mid (R)', 'Frontal Mid Orb (R)', 
# 	'Frontal Sup (R)','Frontal Sup Medial (R)', 'Frontal Sup Orb (R)', 'Hippocampus (R)', 'Insula (R)','ParaHippocampal (R)', 
# 	'Temporal Inf (R)', 'Temporal Mid (R)', 'Temporal Pole Mid (R)','Temporal Pole Sup (R)', 'Temporal Sup (R)'])
# print(len(df2))
sc.add_to_subplot(s_obj_r, row=0, col=1)

# =============================================================================
#						Electrodes sources by subject - LEFT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
b_obj_l.alpha = 0.09
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left',
					title='Electrodes localization by subject - Left')
s_obj_l = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
			alpha=.95, data=data, text=labels)
s_obj_l.set_visible_sources('left')
# df3 = s_obj_l.analyse_sources(roi_obj='aal',distance=15.,keep_only=['Amygdala (L)', 'Cingulum Ant (L)', 
# 	'Cingulum Mid (L)', 'Frontal Inf Orb (L)', 'Frontal Inf Tri (L)', 'Frontal Mid (L)', 
# 	'Frontal Mid Orb (L)', 'Frontal Sup (L)','Frontal Sup Medial (L)', 'Frontal Sup Orb (L)', 
# 	'Hippocampus (L)', 'Insula (L)','ParaHippocampal (L)','Temporal Inf (L)', 'Temporal Mid (L)', 
# 	'Temporal Pole Mid (L)','Temporal Pole Sup (L)', 'Temporal Sup (L)'])
# print(len(df3))
sc.add_to_subplot(s_obj_l, row=0, col=2)

#sc.preview()
sc.screenshot(f_form_save.format(step),autocrop=True,print_size=(6,7))
