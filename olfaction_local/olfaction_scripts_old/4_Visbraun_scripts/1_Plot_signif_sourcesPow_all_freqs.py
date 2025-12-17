from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '01', 'pow', True
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_E_EpiPerf_LowHigh/')
path_mask = join(path_npz, 'masks_visbrain'+th+'/')
path_to_save = join(path_npz, 'Brain_Plots_th'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}_corr.npy') 
f_form_save = join(path_to_save, '{}_proj_min{}_'+feat+'_{}_th{}_corr.png') 
###############################################################################

freqs = ['0_theta','1_alpha','2_beta','3_gamma'] #'3_gamma' '2_theta','3_alpha','4_beta','5_gamma1','6_gamma2'
rois_to_keep = ['OFC']
wins, radius = [1.0], 1.
min_sigs = [1]
methods = ['s_Mai_RL']
cmap = 'jet'
clim = (-100,100)

for win, freq, min_sig, method in product(wins,freqs,min_sigs,methods):
	arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
	#print(arch_sig.files)
	subjects = arch_sig['su_codes']
	power = ((arch_sig['s_elec_pow1']-arch_sig['s_elec_pow0'])/abs(arch_sig['s_elec_pow0']))*100
	perm_thr, da = arch_sig['s_th_'+th], arch_sig['s_da']
	times = arch_sig['s_time'][17:52]-3
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(700, 900)) #largeur, hauteur
	# =============================================================================from visbrain import Brain, Colorbar


freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
wins = [1]
min_sigs = [2]
methods = ['s_Mai_RL'] #'aal_RL','BA'
cmap = 'jet' 
clim = (-100,100)

for min_sig,win, method in product(min_sigs,wins, methods):
	# =============================================================================
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(400, 1000)) #largeur, hauteur
	# ============================================================================

	for i,freq in enumerate(freqs):
		arch_sig = np.load(npz_form.format('All_subjects',freq,'odor'))
		print(arch_sig.files)
		#load data to plot and mask
		mask = np.load(masks_vis_form.format('MAI_RL',min_sig,freq,str(win),th))
		print(mask.dtype, len(mask))
		bad_pow, good_pow = arch_sig['s_elec_pow0'],arch_sig['s_elec_pow1']

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
		#PLOT POWER CHANGE
		data = ((good_pow_max - bad_pow_max)/ bad_pow_max)*100
		#PLOT AUC Score
		#PLOT Time Score
		data = data[np.where(mask==False)]
		#data = idx_all[np.where(mask==False)]
		
		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='bottom',
							title=freq+' - minwin'+str(win)+' th'+th)
		# #Define a source object with a different color by subject
		xyz = arch_sig['s_xyz']
		xyz_sig = xyz[np.where(mask==False)]
		elecs_labels_sig=None
		s_obj_c = SourceObj('modulation', xyz_sig, radius_min=10.,radius_max=11.)
		s_obj_c.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_c, row=i, col=0)

		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=i, col=1, rotate='left')
		s_obj_l = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10.,radius_max=11., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_l.set_visible_sources('left')
		s_obj_l.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_l, row=i, col=1)

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		b_obj_r2 = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10.,radius_max=11.,
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_r2.set_visible_sources('right')
		s_obj_r2.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_r2, row=i, col=2)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(s_obj_r2, cblabel=feat, border=False, **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=3, width_max=200, height_max=600)
sc.screenshot(f_form_save.format(method, min_sig,th),autocrop=True,print_size=(4,10))
