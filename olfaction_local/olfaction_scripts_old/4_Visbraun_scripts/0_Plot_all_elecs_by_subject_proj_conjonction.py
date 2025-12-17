from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product
from np_replace_values_dict import replace_with_dict

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th = '0.01'
freqs = ['0_theta','1_alpha','2_beta','3_gamma']
rois_to_keep = ['ACC','Amg','Amg-PirT','HC','IFG','Ins','MFG','OFC','PHG',
                'SFG','pPirT']
cmap = 'viridis'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/figure/'
path_npz_enc = join(path, '0_Classif_Power_E_EpiPerf_LowHigh_1000perm_BBG/')
path_npz_ret = join(path, '0_Classif_Power_R_EpiPerf_LowHigh_1000perm_BBG/')
npz_form = 'All_subjects_sources_{}_odor_low_high_sel_physFT.npz'
mask_form = 'All_subjects_mask_stat_{}_minwin1.0_th'+th+'.npy'
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
	
	#select all electrodes from encoding (specific regions and subjects)
	mat_enc = np.load(path_npz_enc+npz_form.format(freq))
	idx_E = np.where([roi in rois_to_keep for roi in mat_enc['s_labels']])[0]
	id_su = np.where([su != 'S6' for su in mat_enc['su_codes'][idx_E]])[0]
	xyz0, labels0 = mat_enc['s_xyz'][idx_E][id_su], mat_enc['s_labels'][idx_E][id_su]
	mask0 = np.load(path_npz_enc+'masks_stat/'+mask_form.format(freq))[id_su]
	mask0 = np.asarray([not x for x in mask0])
	    
	#select all electrodes from retrieval (specific regions, all subjects)
	mat_ret = np.load(path_npz_ret+npz_form.format(freq))
	idx_R = np.where([roi in rois_to_keep for roi in mat_ret['s_labels']])[0]
	xyz1, labels1 = mat_ret['s_xyz'][idx_R], mat_ret['s_labels'][idx_R]
	mask1 = np.load(path_npz_ret+'masks_stat/'+mask_form.format(freq))
	mask1 = np.asarray([not x for x in mask1])

	#select only electrodes present both during encoding and retrieval
	to_keep_enc = np.where((xyz0==xyz1[:,None]).all(-1))[1]
	to_keep_ret = np.where((xyz1==xyz0[:,None]).all(-1))[1]
	mask0_new, mask1_new = mask0[to_keep_enc], mask1[to_keep_ret]
	print(len(mask0_new),len(mask1_new))

	xyz_shared = xyz0[to_keep_enc]
	xyz0_sig = xyz0[to_keep_enc][mask0_new]
	xyz1_sig = xyz1[to_keep_ret][mask1_new]
	print(freq,xyz0_sig.shape,xyz1_sig.shape)

	#loop on all elecs and create a code indicating when it decodes
	codes_elecs = []
	for i, elec in enumerate(xyz_shared):
		in_1 = np.where(np.all(elec == xyz0_sig, 1))[0] #encoding
		in_2 = np.where(np.all(elec == xyz1_sig, 1))[0] #retrieval

		if in_1.size and in_2.size: #decode in both
		    codes_elecs.append(2)
		if in_1.size and not in_2.size: #decode only in encoding
		    codes_elecs.append(0)
		if not in_1.size and in_2.size: #decode only in retrieval
		    codes_elecs.append(1)
		if not in_1.size and not in_2.size:
			codes_elecs.append(-1)

	print(codes_elecs)
	#loop on all elecs and create a code indicating when it decodes
	# -1 when elec is not signif whatever the condition
	# 0 when elec is signif only during encoding
	# 1 when elec is signif only during retrieval
	# 2 when elec is signif during encoding and retrieval

	# =============================================================================
	#						Electrodes sources by subject - RIGHT EXT 
	# =============================================================================
	u_color = {-1:"grey", 0:"#114693", 1:"#30169A",2:"#DA0615"} #encoding bleu, retrieval orange, conj jaune
	mask_elec = [False if code == -1 else True for code in codes_elecs]
	codes_elecs = np.array(codes_elecs)[mask_elec]
	color = [u_color[k] for k in codes_elecs]
	xyz_all = xyz_shared[mask_elec]
	data = np.arange(len(codes_elecs))

	#Define a brain object to plot
	b_obj_0 = BrainObj('B2', translucent=False)
	b_obj_0.alpha = 0.09
	sc.add_to_subplot(b_obj_0, row=0, col=0, rotate='right')
	s_obj_c = SourceObj('Color', xyz_all, symbol='disc', color=color,
		radius_min=10., alpha=.95, data=codes_elecs)
	b_obj_0.project_sources(s_obj_c,'modulation', cmap=cmap,radius=10.)
	s_obj_c.visible_obj = False
	sc.add_to_subplot(s_obj_c, row=0, col=0)
	#sc.preview()

	# =============================================================================
	#						Electrodes sources by subject - RIGHT INT 
	# =============================================================================
	#Define a brain object to plot
	b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
	b_obj_r.alpha = 0.09
	sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='left')
	s_obj_r = SourceObj('S_right', xyz_all, symbol='disc', color=color,radius_min=10., 
				alpha=.95, data=codes_elecs)
	b_obj_r.project_sources(s_obj_r,'modulation',cmap=cmap,radius=10.)
	s_obj_r.visible_obj = False
	sc.add_to_subplot(s_obj_r, row=0, col=1)

	# =============================================================================
	#						Electrodes sources by subject - LEFT EXT 
	# =============================================================================
	#Define a brain object to plot
	b_obj_l = BrainObj('B2', translucent=False)
	b_obj_l.alpha = 0.09
	sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left')
	s_obj_l = SourceObj('S_right', xyz_all, symbol='disc', color=color,radius_min=10., 
				alpha=.95, data=codes_elecs)
	b_obj_l.project_sources(s_obj_l,'modulation',cmap=cmap,radius=10.)
	s_obj_l.visible_obj = False
	sc.add_to_subplot(s_obj_l, row=0, col=2)

	# =============================================================================
	#						Electrodes sources by subject - LEFT INT 
	# =============================================================================
	#Define a brain object to plot
	b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='left')
	b_obj_l2.alpha = 0.09
	sc.add_to_subplot(b_obj_l2, row=0, col=3, rotate='right')
	s_obj_l2 = SourceObj('S_right', xyz_all, symbol='disc', color=color,radius_min=10., 
				alpha=.95, data=codes_elecs)
	b_obj_l2.project_sources(s_obj_l2,'modulation',cmap=cmap,radius=10.)
	s_obj_l2.visible_obj = False
	sc.add_to_subplot(s_obj_l2, row=0, col=3)

	cb_rep = ColorbarObj(b_obj_l2, cblabel='Decoding', **CBAR_STATE)
	sc.add_to_subplot(cb_rep, row=0, col=4, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(freq),autocrop=True,print_size=(6,7))
