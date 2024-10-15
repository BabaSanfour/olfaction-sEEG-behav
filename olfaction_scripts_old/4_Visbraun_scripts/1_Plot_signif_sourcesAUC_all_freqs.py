from visbrain import Brain, Colorbar
from os import listdir, makedirs
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '05', 'powAUC', True
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_E_EpiPerf_LowHigh/')
path_mask = join(path_npz, 'masks_visbrain'+th+'/')
path_to_save = join(path_npz, 'Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}_corr.npy')
f_form_save = join(path_to_save, '{}_sources_min{}_'+feat+'_allfreqs_th{}_corr.png')
###############################################################################

freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
wins = [1]
min_sigs = [1]
methods = ['s_Mai_RL'] #'aal_RL','BA'
cmap = 'autumn' 
clim = (0.7,1)

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

		#power change corresponding to max AUC score
		da = arch_sig['s_da']
		perm_max = round(max(arch_sig['s_th_'+th]),2)
		da_max = np.array([])
		for elec in range(da.shape[0]):
			da_elec = max(da[elec])
			da_max = np.hstack((da_max, da_elec)) if np.size(da_max) else da_elec
		data = da_max
		data = data[np.where(mask==False)]
		
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
		cb_rep = ColorbarObj(s_obj_r2, cblabel=feat, border=False, **CBAR_STATE,vmin=perm_max, under='grey')
		sc.add_to_subplot(cb_rep, row=i, col=3, width_max=200, height_max=600)
sc.screenshot(f_form_save.format(method, min_sig,th),autocrop=True,print_size=(4,10))
