from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj,RoiObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '0.01', 'pow', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_E_Odor_Inspi_1000perm/')
path_mask = join(path_npz, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
f_form_save = join(path_to_save, '{}_proj_LTM_min{}_win{}_'+feat+'_{}_th{}.png') 
###############################################################################

freqs = ['2_theta','3_alpha','4_beta','5_gamma1','6_gamma2']
list_su = ['S0','S1','S2','S3','S4','S5','S6']
wins, radius = [1], 10.
min_sigs = [1]
methods = ['s_Mai_RL']
clim = (-100,100)

for win, min_sig, method in product(wins,min_sigs, methods):
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(300,450)) #largeur, hauteur
	# ============================================================================

	for i,freq in enumerate(freqs):
		arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
		da, power = arch_sig['s_da'], ((arch_sig['s_elec_pow1']-arch_sig['s_elec_pow0'])/abs(arch_sig['s_elec_pow0']))*100
		if exists(masks_form.format(freq,str(win),th)):
			mask = np.load(masks_form.format(freq,str(win),th))
		else:
			mask = np.ones(da.shape[0])#create a mask for all elecs when non significant
		xyz = arch_sig['s_xyz']
		#Find the max AUC score and corresponding power
		idx_all, pow_all = np.array([]), np.array([])
		for elec in range(da.shape[0]):
			da_elec = da[elec]
			idx = [i for i,j in enumerate(da_elec) if j ==max(da_elec)][0]
			pow = power[elec,idx]
			pow_all = np.hstack((pow_all,pow)) if np.size(pow_all) else pow
			idx_all = np.hstack((idx_all, idx)) if np.size(idx_all) else idx
		data = pow_all
		print(data.shape, idx_all.shape)

		#Select only electrodes in the MTL
		rois = arch_sig['s_MAI_RL']
		sel_idx = []
		for e in range(rois.shape[0]):
			if rois[e] in ['HC','Amg','PHG','pPirT','Amg-PirT']:
				sel_idx.append(e)
		xyz = xyz[sel_idx]
		data = data[sel_idx]
		mask = mask[sel_idx]
				
		# =============================================================================
		#						Electrodes sources by subject - Front VIEW 
		# =============================================================================
		b_obj_f = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_f.alpha = 0.09
		sc.add_to_subplot(b_obj_f, row=i, col=0, rotate='front')
		s_obj_f = SourceObj('S_top', xyz, data=data, mask=mask)
		s_obj_f.visible_obj = False
		sc.add_to_subplot(s_obj_f, row=i, col=0)
		roi_obj_f = RoiObj('aal', border=False)
		idx = roi_obj_f.where_is(['Hippocampus','Amygdala','ParaHippocampal'])
		roi_obj_f.select_roi(select=idx, unique_color=True, smooth=10)
		roi_obj_f.project_sources (s_obj_f, 'modulation', cmap='jet', clim=clim, radius=radius)
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
		roi_obj = RoiObj('aal', border=False)
		idx = roi_obj.where_is(['Hippocampus','Amygdala','ParaHippocampal'])
		roi_obj.select_roi(select=idx, unique_color=True, smooth=10)
		roi_obj.project_sources (s_obj_r, 'modulation', cmap='jet', clim=clim, radius=radius)
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
		roi_obj_l = RoiObj('aal', border=False)
		idx = roi_obj_l.where_is(['Hippocampus','Amygdala','ParaHippocampal'])
		roi_obj_l.select_roi(select=idx, unique_color=True, smooth=10)
		roi_obj_l.project_sources (s_obj_l, 'modulation', cmap='jet', clim=clim, radius=radius)
		s_obj_l.fit_to_vertices(roi_obj_l.vertices)
		sc.add_to_subplot(roi_obj_l, row=i, col=2)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(roi_obj_l, cblabel='ΔPower (%)', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=3, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method,min_sig,win,'allfreqs',th),autocrop=True,print_size=(9,12))