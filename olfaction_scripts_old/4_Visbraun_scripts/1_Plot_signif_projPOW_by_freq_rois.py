from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj, RoiObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '0.01', 'pow', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/LDA_TPSim_E_R_all_BBG/')
path_mask = join(path_npz, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
f_form_save = join(path_to_save, '{}_proj_min{}_win{}_'+feat+'_{}_th{}_{}.png') 
###############################################################################

freqs = ['0_theta','1_alpha','2_beta','3_gamma'] #'3_gamma' '2_theta','3_alpha','4_beta','5_gamma1','6_gamma2'
rois_to_keep = ['OFC']
wins, radius = [1.0], 1.
min_sigs = [1]
methods = ['s_Mai_RL']

for win, min_sig, method in product(wins,min_sigs, methods):
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(300,300)) #largeur, hauteur
	# ============================================================================

	#for i,freq in zip(range(2,7),freqs):
	for i,freq in enumerate(freqs):
		clim = (0,100)

		arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
		#print(arch_sig.files,arch_sig['s_labels'][:10],arch_sig['s_MAI_RL'][:10])
		idx_rois = np.where([roi in rois_to_keep for roi in arch_sig['s_labels']])[0]
		da = arch_sig['s_da'][idx_rois]#,10:30]#[idx_rois,:][0]
		print('da shape', da.shape)

		subjects = arch_sig['su_codes'][idx_rois]
		xyz = arch_sig['s_xyz'][idx_rois]
		pow0, pow1 = arch_sig['s_elec_pow0'][idx_rois], arch_sig['s_elec_pow1'][idx_rois]
		print('pow shape',pow0.shape,pow1.shape)

		power = ((pow1-pow0)/abs(pow0))*100
		mask = np.load(masks_form.format(freq,str(win),th))[idx_rois]
		mask = np.asarray([x if power[i] > 0  else True for i,x in enumerate(mask)])
		data = power
		print(data.shape,xyz.shape,mask.shape)#idx_all.shape)

				
		# =============================================================================
		#						Electrodes sources by subject - Front VIEW 
		# =============================================================================
		b_obj_f = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_f.alpha = 0.09
		sc.add_to_subplot(b_obj_f, row=i, col=0, rotate='front')
		s_obj_f = SourceObj('S_top', xyz, data=data,mask=mask)
		s_obj_f.visible_obj = False
		sc.add_to_subplot(s_obj_f, row=i, col=0)
		roi_obj_f = RoiObj('aal',cblabel='TPSim High/Low', border=False)
		idx = roi_obj_f.where_is(['Hippocampus'])
		roi_obj_f.select_roi(select=idx, smooth=None)
		roi_obj_f.project_sources (s_obj_f, 'modulation', cmap='jet', clim=clim, radius=radius,
			under='gray')
		s_obj_f.fit_to_vertices(roi_obj_f.vertices)
		sc.add_to_subplot(roi_obj_f, row=i, col=0)

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=1, rotate='right',title=freq+' - win'+str(win)+' th'+th)
		s_obj_r = SourceObj('S_right', xyz, data=data, mask=mask)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=1)
		roi_obj = RoiObj('aal', cblabel='TPSim High/Low',border=False)
		idx = roi_obj.where_is(['Hippocampus'])
		roi_obj.select_roi(select=idx, smooth=3.)
		roi_obj.project_sources (s_obj_r, 'modulation', cmap='jet', clim=clim, radius=radius,
			under='gray')
		s_obj_r.fit_to_vertices(roi_obj.vertices)
		sc.add_to_subplot(roi_obj, row=i, col=1)

		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=i, col=2, rotate='left')
		s_obj_l = SourceObj('S_left', xyz, data=data, mask=mask)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=2)
		roi_obj_l = RoiObj('aal', cblabel='TPSim High/Low',border=False)
		idx = roi_obj_l.where_is(['Hippocampus'])
		roi_obj_l.select_roi(select=idx, smooth=None)
		roi_obj_l.project_sources (s_obj_l, 'modulation', cmap='jet', clim=clim, radius=radius,
			under='gray')
		s_obj_l.fit_to_vertices(roi_obj_l.vertices)
		sc.add_to_subplot(roi_obj_l, row=i, col=2)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(roi_obj_l, cblabel='ΔPower (%)', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=3, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method,min_sig,win,'allfreqs',th,rois_to_keep[0]),
						autocrop=True,print_size=(9,12))