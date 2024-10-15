from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product
from np_replace_values_dict import replace_with_dict

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th = '0.01'
freqs = ['2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/figure/'
path_npz_enc = join(path, '0_Classif_Power_E_EpiPerf_LowHigh_1000perm/')
path_npz_ret = join(path, '0_Classif_Power_R_EpiPerf_LowHigh_1000perm/')
npz_form = 'All_subjects_sources_{}_odor_low_high_sel_physFT.npz'
mask_form = 'All_subjects_mask_stat_{}_minwin1_th'+th+'.npy'
f_form_save = join(path, 'Conjonction/Conj_distribution_proj_{}.png')
###############################################################################


for freq in freqs:
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(700, 500)) #largeur, hauteur
	# =============================================================================
	#					Electrodes codes and xyz - To Plot
	# =============================================================================

	mat_enc = np.load(path_npz_enc+npz_form.format(freq))
	xyz0 = mat_enc['s_xyz']
	mask0 = np.load(path_npz_enc+'masks_stat/'+mask_form.format(freq))
	mask_sig0 = [False if i == True else True for i in mask0]	
	mat_ret = np.load(path_npz_ret+npz_form.format(freq))
	xyz1 = mat_ret['s_xyz']
	mask1 = np.load(path_npz_ret+'masks_stat/'+mask_form.format(freq))
	mask_sig1 = [False if i == True else True for i in mask1]
	
	#Select only significant electrodes
	xyz_sig0 = xyz0[mask_sig0]
	xyz_sig1 = xyz1[mask_sig1]
	data0 = np.ones(xyz_sig0.shape[0])
	data1 = np.ones(xyz_sig1.shape[0])

	# =============================================================================
	#						Electrodes sources by subject - TOP VIEW 
	# =============================================================================
	#Define a brain object to plot
	b_obj_top = BrainObj('B2', translucent=False)
	b_obj_top.alpha = 0.09
	sc.add_to_subplot(b_obj_top, row=0, col=0)
	s_obj_c0 = SourceObj('Color', xyz_sig0, data=data0)
	s_obj_c1 = SourceObj('Color', xyz_sig1, data=data1)
	b_obj_top.project_sources(s_obj_c0, cmap='RdPu', to_overlay=0, clim=(0., 1.))
	b_obj_top.project_sources(s_obj_c1, cmap='Wistia', to_overlay=1, clim=(0., 1.))
	s_obj_c0.visible_obj = False
	s_obj_c1.visible_obj = False
	sc.add_to_subplot(s_obj_c0, row=0, col=0)
	sc.add_to_subplot(s_obj_c1, row=0, col=0)
	#sc.preview()

	# =============================================================================
	#						Electrodes sources by subject - RIGHT VIEW 
	# =============================================================================
	#Define a brain object to plot
	b_obj_r = BrainObj('B2', translucent=False)
	b_obj_r.alpha = 0.09
	sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
	s_obj_r0 = SourceObj('S_right', xyz_sig0, data=data0)
	s_obj_r1 = SourceObj('S_right', xyz_sig1, data=data1)
	b_obj_r.project_sources(s_obj_r0,cmap='RdPu', to_overlay=0, clim=(0., 1.))
	b_obj_r.project_sources(s_obj_r1,cmap='Wistia', to_overlay=1, clim=(0., 1.))
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
	b_obj_l.project_sources(s_obj_l0,cmap='RdPu', to_overlay=0, clim=(0., 1.))
	b_obj_l.project_sources(s_obj_l1,cmap='Wistia', to_overlay=1, clim=(0., 1.))
	s_obj_l0.visible_obj = False
	s_obj_l1.visible_obj = False
	sc.add_to_subplot(s_obj_l0, row=0, col=2)
	sc.add_to_subplot(s_obj_l1, row=0, col=2)

	#sc.preview()
	sc.screenshot(f_form_save.format(freq),autocrop=True,print_size=(6,7))
