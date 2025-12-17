from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

score, feat = 'Epi', 'POS'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_odor_expi/')
path_mask = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score/masks_visbrain/')
path_to_save = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score/Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_'+score+'Score_Expi_AAL_BA_{}.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_'+feat+'_bsl_{}_{}_3wins_th{}.png')
###############################################################################

freqs = ['2_theta', '3_alpha','4_beta','5_gamma1','6_gamma2'] #'0_VLFC', '1_delta', 
bsls = ['None']
wins = [3]
ths = ['05']
min_sigs = [4] #for min=0 all methods are equivalent
methods = ['aal_RL']#'aal','aal_RL',BA

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
		mask = np.load(masks_vis_form.format(method,min_sig,freq,str(win),th))
		#print(mask)

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='top',
							title=freq+' win '+str(win)+' th'+th+' - top')
		#Define a source object with a different color by subject
		xyz, subjects, elecs, channels = arch_sig['s_xyz'],arch_sig['su_codes'], arch_sig['s_elec'], arch_sig['s_channels']
		xyz_sig = xyz[np.where(mask==False)]
		subjects_sig = subjects[np.where(mask==False)]
		#elecs_sig, channels_sig = elecs[np.where(mask==False)], channels[np.where(mask==False)]
		#elecs_labels_sig = [elecs_sig[i] for i in range(xyz_sig.shape[0])]
		elecs_labels_sig=None
		u_color = ["darkblue", "royalblue", "deepskyblue", "mediumspringgreen", "yellow", "darkorange","red"]
		color = [u_color[int(k[1])] for k in subjects_sig]
		data = np.arange(len(subjects_sig))
		s_obj_c = SourceObj('Color', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		sc.add_to_subplot(s_obj_c, row=i, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=1, rotate='left',
							title=freq+' win '+str(win)+' th'+th+' - Right')
		s_obj_r = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_r.set_visible_sources('right')
		sc.add_to_subplot(s_obj_r, row=i, col=1)

		b_obj_r2 = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_r2.set_visible_sources('right')
		sc.add_to_subplot(s_obj_r2, row=i, col=2)

		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=i, col=3, rotate='left',
							title=freq+' win '+str(win)+' th'+th+' - Left')
		s_obj_l = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_l.set_visible_sources('left')
		sc.add_to_subplot(s_obj_l, row=i, col=3)

		b_obj_l2 = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_l2.set_visible_sources('left')
		sc.add_to_subplot(s_obj_l2, row=i, col=4)

	#sc.preview()
	sc.screenshot(f_form_save.format(method, min_sig,bsl,freq,th),autocrop=True,print_size=(9,12))