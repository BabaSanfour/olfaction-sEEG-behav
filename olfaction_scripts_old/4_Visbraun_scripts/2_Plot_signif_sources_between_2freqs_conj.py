from visbrain import Brain, Colorbar
from os import listdir, makedirs
from os.path import isfile, join, exists
import numpy as np
from itertools import product
from np_replace_values_dict import replace_with_dict

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th = '0.01'
freqs = ['2_theta','3_alpha', '6_gamma2'] #'3_alpha', '4_beta','5_gamma1',
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/figure/'
path_npz_enc = join(path, '0_Classif_Power_E_EpiPerf_LowHigh_1000perm/')
npz_form = 'All_subjects_sources_{}_odor_low_high_sel_physFT.npz'
mask_form = 'All_subjects_mask_stat_{}_minwin1_th'+th+'.npy'
path2save = join(path_npz_enc, 'Freqs_conjonction/')
f_form_save = join(path2save, 'Conj_distribution_proj_{}.png')
###############################################################################
if not exists(path2save):
	makedirs(path2save)
###############################################################################


xyz_sig_freqs, data_freqs = [],[]
for freq in freqs:

	#Load all electrodes coordinates
	mat = np.load(path_npz_enc+npz_form.format(freq))
	xyz = mat['s_xyz']
	mask = np.load(path_npz_enc+'masks_stat/'+mask_form.format(freq))
	mask_sig = [False if i == True else True for i in mask]	
		
	#Select only significant electrodes
	xyz_sig = xyz[mask_sig]
	data = np.ones(xyz_sig.shape[0])
	xyz_sig_freqs.append(xyz_sig)
	data_freqs.append(data)
print([xyz.shape for xyz in xyz_sig_freqs], [data.shape for data in data_freqs])

# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(700, 500)) #largeur, hauteur

# =============================================================================
#						Electrodes sources by subject - TOP VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_top = BrainObj('B2', translucent=False)
b_obj_top.alpha = 0.09
sc.add_to_subplot(b_obj_top, row=0, col=0)
s_obj_c0 = SourceObj('Color', xyz_sig_freqs[0], data=data_freqs[0])
s_obj_c1 = SourceObj('Color', xyz_sig_freqs[1], data=data_freqs[1])
s_obj_c2 = SourceObj('Color', xyz_sig_freqs[2], data=data_freqs[2])
b_obj_top.project_sources(s_obj_c0, cmap='Reds', to_overlay=0, clim=(0., 1.))
b_obj_top.project_sources(s_obj_c1, cmap='Blues', to_overlay=1, clim=(0., 1.))
b_obj_top.project_sources(s_obj_c2, cmap='spring', to_overlay=2, clim=(0., 1.))
s_obj_c0.visible_obj = False
s_obj_c1.visible_obj = False
s_obj_c2.visible_obj = False
sc.add_to_subplot(s_obj_c0, row=0, col=0)
sc.add_to_subplot(s_obj_c1, row=0, col=0)
sc.add_to_subplot(s_obj_c2, row=0, col=0)
sc.preview()
sc.screenshot(f_form_save.format(freq),autocrop=True,print_size=(6,7))


# =============================================================================
#						Electrodes sources by subject - RIGHT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=False)
b_obj_r.alpha = 0.09
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
s_obj_r0 = SourceObj('S_right', xyz_sig0, data=data0)
s_obj_r1 = SourceObj('S_right', xyz_sig1, data=data1)
b_obj_r.project_sources(s_obj_r0,cmap='Reds', to_overlay=0, clim=(0., 1.))
b_obj_r.project_sources(s_obj_r1,cmap='Blues', to_overlay=1, clim=(0., 1.))
s_obj_r0.visible_obj = False
s_obj_r1.visible_obj = False
sc.add_to_subplot(s_obj_r0, row=0, col=1)
sc.add_to_subplot(s_obj_r1, row=0, col=1)

# =============================================================================
#						Electrodes sources by subject - LEFT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=False)
b_obj_l.alpha = 0.09
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
s_obj_l0 = SourceObj('S_right', xyz_sig0, data=data0)
s_obj_l1 = SourceObj('S_right', xyz_sig1, data=data1)
b_obj_l.project_sources(s_obj_l0,cmap='Reds', to_overlay=0, clim=(0., 1.))
b_obj_l.project_sources(s_obj_l1,cmap='Blues', to_overlay=1, clim=(0., 1.))
s_obj_l0.visible_obj = False
s_obj_l1.visible_obj = False
sc.add_to_subplot(s_obj_l0, row=0, col=2)
sc.add_to_subplot(s_obj_l1, row=0, col=2)

#sc.preview()
sc.screenshot(f_form_save.format(freq),autocrop=True,print_size=(6,7))
