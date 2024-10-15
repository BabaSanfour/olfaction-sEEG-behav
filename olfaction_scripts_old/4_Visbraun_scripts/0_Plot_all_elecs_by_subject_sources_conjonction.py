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
rois_to_keep = ['ACC','Amg','Amg-PirT','HC','IFG','Ins','MFG','OFC','PHG',
                'SFG','pPirT']
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/figure/'
path_npz_enc = join(path, '0_Classif_Power_E_EpiPerf_LowHigh_1000perm_v2/')
path_npz_ret = join(path, '0_Classif_Power_R_EpiPerf_LowHigh_1000permv2/')
npz_form = 'All_subjects_sources_{}_odor_low_high_sel_physFT.npz'
mask_form = 'All_subjects_mask_stat_{}_minwin1.0_th'+th+'.npy'
f_form_save = join(path, 'Conjonction/Conj_distribution_{}.png')
###############################################################################
radius_min = 10.

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
	idx_E = np.where([roi in rois_to_keep for roi in mat_enc['s_MAI_RL']])
	xyz0 = mat_enc['s_xyz'][idx_E]
	mask0 = np.load(path_npz_enc+'masks_stat/'+mask_form.format(freq))

	mat_ret = np.load(path_npz_ret+npz_form.format(freq))
	idx_R = np.where([roi in rois_to_keep for roi in mat_ret['s_MAI_RL']])
	xyz1 = mat_ret['s_xyz'][idx_R]
	mask1 = np.load(path_npz_ret+'masks_stat/'+mask_form.format(freq))

	xyz_all = np.unique(np.concatenate((xyz0,xyz1),axis=0),axis=0)
	print(xyz0.shape,xyz1.shape,xyz_all.shape)

	#loop on all elecs and create a code indicating when it decodes
	# -1 when elec is not signif whatever the condition
	# 0 when elec is signif only during encoding
	# 1 when elec is signif only during retrieval
	# 2 when elec is signif during encoding and retrieval

	codes_elecs = []
	for i, elec in enumerate(xyz_all):
		in_1 = np.where(np.all(elec == xyz0, 1))[0] #encoding
		in_2 = np.where(np.all(elec == xyz1, 1))[0] #retrieval

		if in_1.size and not in_2.size:
			sig = mask0[in_1]
			codes_elecs.append([-1 if sig == True else 0][0]) #True means masked so NOT signif
		if not in_1.size and in_2.size:
			sig = mask1[in_2]
			codes_elecs.append([-1 if sig == True else 1][0])
		if in_1.size and in_2.size:
			sig0,sig1 = mask0[in_1], mask1[in_2]
			if sig0 == False and sig1 == False:
				codes_elecs.append(2)
			if sig0 == True and sig1 == True:
				codes_elecs.append(-1)
			if sig0 == False and sig1 == True:
				codes_elecs.append(0)
			if sig0 == True and sig1 == False:
				codes_elecs.append(1)

	# =============================================================================
	#						Electrodes sources by subject - TOP VIEW 
	# =============================================================================
	#u_color = {-1:"grey", 0:"#114693", 1:"#DA0615",2:"#30169A"} #0:"#0B2FB2", 1:"#FF7F00",2:"#FFDA00"#encoding bleu, retrieval orange, conj jaune
	u_color = {-1:"grey", 0:"indigo", 1:"yellow",2:"teal"} 
	mask_elec = [False if code == -1 else True for code in codes_elecs]
	codes_elecs = np.array(codes_elecs)[mask_elec]
	color = [u_color[k] for k in codes_elecs]
	xyz_all = xyz_all[mask_elec]
	data = np.arange(len(codes_elecs))

	#Define a brain object to plot
	b_obj_top = BrainObj('B2', translucent=True)
	b_obj_top.alpha = 0.09
	sc.add_to_subplot(b_obj_top, row=0, col=0)
	s_obj_c = SourceObj ('Color', xyz_all, symbol='disc', color=color,
		radius_min=radius_min, alpha=.95, data=data)
	sc.add_to_subplot(s_obj_c, row=0, col=0)
	#sc.preview()

	# =============================================================================
	#						Electrodes sources by subject - RIGHT VIEW 
	# =============================================================================
	#Define a brain object to plot
	b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
	b_obj_r.alpha = 0.09
	sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right')
	s_obj_r = SourceObj('S_right', xyz_all, symbol='disc', color=color,radius_min=radius_min, 
				alpha=.95, data=data)
	s_obj_r.set_visible_sources('right')
	sc.add_to_subplot(s_obj_r, row=0, col=1)

	# =============================================================================
	#						Electrodes sources by subject - LEFT VIEW 
	# =============================================================================
	#Define a brain object to plot
	b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
	b_obj_l.alpha = 0.09
	sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
	s_obj_l = SourceObj('S_right', xyz_all, symbol='disc', color=color,radius_min=radius_min, 
				alpha=.95, data=data)
	s_obj_l.set_visible_sources('left')
	sc.add_to_subplot(s_obj_l, row=0, col=2)

	#sc.preview()
	sc.screenshot(f_form_save.format(freq),autocrop=True,print_size=(6,7))
