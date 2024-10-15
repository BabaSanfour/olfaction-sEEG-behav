from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj,RoiObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '0.001', 'pow', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/LDA_Power_R_pre_stim_1s/')
path_mask = join(path_npz, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
f_form_save = join(path_to_save, '{}_proj_min{}_win{}_'+feat+'_{}_th{}.png') 
###############################################################################

freqs = ['1_delta','2_theta','3_alpha','4_beta','5_gamma1','6_gamma2'] #'0_VLFC' (none at 0.001)
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
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(700, 900)) #largeur, hauteur
	# ============================================================================

	#for i,freq in zip(range(2,7),freqs):
	for i,freq in enumerate(freqs):
		clim = (-80,80) if freq in ['4_beta','5_gamma1','6_gamma2'] else (-100,100) 
		arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
		da = arch_sig['s_da']
		power = ((arch_sig['s_elec_pow1']-arch_sig['s_elec_pow0'])/abs(arch_sig['s_elec_pow0']))*100
		subjects = arch_sig['su_codes']
		if exists(masks_form.format(freq,str(win),th)):
			mask = np.load(masks_form.format(freq,str(win),th))
		else:
			mask = np.ones(da.shape[0])#create a mask for all elecs when non significant
		
		#Find the max AUC score and corresponding power
		idx_all, pow_all = np.array([]), np.array([])
		for elec in range(da.shape[0]):
			da_elec = da[elec]
			idx = [i for i,j in enumerate(da_elec) if j ==max(da_elec)][0]
			pow = power[elec,idx]
			pow_all = np.hstack((pow_all,pow)) if np.size(pow_all) else pow
			idx_all = np.hstack((idx_all, idx)) if np.size(idx_all) else idx
		data = pow_all #power pow_all
		#print(data.shape, idx_all.shape)

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='bottom',
							title=freq+' - win'+str(win)+' th'+th, title_color='black')
		#Define a source object with a different color by subject
		xyz = arch_sig['s_xyz']
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_c.visible_obj = False
		sc.add_to_subplot(s_obj_c, row=i, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=1, rotate='left')
		s_obj_r = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r.project_sources(s_obj_r,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=1)

		b_obj_r2 = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_r2.visible_obj = False
		sc.add_to_subplot(s_obj_r2, row=i, col=2)
		
		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=i, col=3, rotate='left')
		s_obj_l = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l.project_sources(s_obj_l,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=3)

		b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=i, col=4)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='ΔTPSim (%)', txtcolor='black', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=5, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method,min_sig,win,'allfreqs',th),autocrop=True,print_size=(9,12))