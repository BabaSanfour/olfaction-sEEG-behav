from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

score, feat = 'Epi', 'TIME'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_odor_expi_AAL_selection/')
path_mask = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score_4aal/masks_visbrain/')
path_to_save = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score_4aal/Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_'+score+'Score_Expi_sel_aal_4_{}.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_'+feat+'_bsl_{}_{}_3wins_th{}.png')
###############################################################################

freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha','4_beta','5_gamma1','6_gamma2'] #
bsls = ['None']
wins = [3,4,5]
ths = ['05']
min_sigs = [0] #for min=0 all methods are equivalent
methods = ['aal_RL']#'aal','aal_RL',BA
clim_min, clim_max = 0, 43
cmap = 'winter'

for freq, bsl, th, min_sig, method in product(freqs,bsls,ths,min_sigs, methods):
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(1000, 700)) #largeur, hauteur
	# ============================================================================
	arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor','None'))
	#print(arch_sig.files)
	
	for i,win in enumerate(wins):
		#load data to plot and mask
		mask = np.load(masks_vis_form.format(method,min_sig,freq,str(win),th))
		bad_pow, good_pow = arch_sig['s_elec_pow_bad'],arch_sig['s_elec_pow_good']

		#compute power change to plot
		#mean power change for significant electrodes across time 4.5s
		#data = ((np.mean(good_pow, axis=1)-np.mean(bad_pow, axis=1))/np.mean(bad_pow, axis=1))*100
		
		#power change corresponding to max AUC score
		da = arch_sig['s_da']
		bad_pow_max, good_pow_max, idx_all = np.array([]), np.array([]), np.array([])
		for elec in range(da.shape[0]):
			da_elec = da[elec]
			idx = [i for i,j in enumerate(da_elec) if j ==max(da_elec)]
			bad_pow_elec, good_pow_elec = np.mean(bad_pow[elec][idx]), np.mean(good_pow[elec][idx])
			#print(bad_pow_elec, good_pow_elec)
			bad_pow_max = np.hstack((bad_pow_max,bad_pow_elec)) if np.size(bad_pow_max) else bad_pow_elec
			good_pow_max = np.hstack((good_pow_max, good_pow_elec)) if np.size(good_pow_max) else good_pow_elec
			idx_all = np.hstack((idx_all, np.mean(idx))) if np.size(idx_all) else idx
		#data = ((good_pow_max - bad_pow_max) / bad_pow_max)*100
		data = idx_all

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B3', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='top',
							title=freq+' - minwin'+str(win)+' th'+th)
		#Define a source object with a different color by subject
		xyz = arch_sig['s_xyz']
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_c.visible_obj = False
		sc.add_to_subplot(s_obj_c, row=i, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - BOTTOM VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_bot = BrainObj('B3', translucent=False, hemisphere='both')
		b_obj_bot.alpha = 0.09
		sc.add_to_subplot(b_obj_bot, row=i, col=1, rotate='bottom')
		s_obj_b = SourceObj('bottom', xyz, mask=mask,data=data)
		b_obj_bot.project_sources(s_obj_b,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_b.visible_obj = False
		sc.add_to_subplot(s_obj_b, row=i, col=1)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B3', translucent=False, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=2, rotate='left',title='Right hemisphere')
		s_obj_r = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r.project_sources(s_obj_r,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=2)

		b_obj_r2 = BrainObj('B3', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=3, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_r2.visible_obj = False
		sc.add_to_subplot(s_obj_r2, row=i, col=3)
		
		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B3', translucent=False, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=i, col=4, rotate='left', title='Left hemisphere')
		s_obj_l = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l.project_sources(s_obj_l,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=4)

		b_obj_l2 = BrainObj('B3', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=5, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap=cmap, clim=(clim_min,clim_max),radius=12.)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=i, col=5)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel=feat, border=False,**CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=6, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method, min_sig,bsl,freq,th),autocrop=True,print_size=(9,12))