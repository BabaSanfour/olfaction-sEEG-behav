from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_odor_expi/Power//')
path_mask = join(path, 'Visbrain/All_subjects/allrois_time/masks_visbrain/')
path_to_save = join(path, 'Visbrain/All_subjects/allrois_time/Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_EpiScore_Expi_AAL_BA_{}.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_POW_bsl_{}_{}_3wins_th{}.png')
###############################################################################

freqs = ['5_gamma1','6_gamma2']#'0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta'
bsls = ['None']
wins = [2,3,4]
ths = ['05']
min_sigs = [3,4]
methods = ['aal','aal_RL','BA']#'aal','aal_RL',BA

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
	
	for win in wins:
		#load data to plot and mask
		mask = np.load(masks_vis_form.format(method,min_sig,freq,str(win),th))
		data = ((arch_sig['s_elec_pow_good'][:,win]-arch_sig['s_elec_pow_bad'][:,win])/arch_sig['s_elec_pow_bad'][:,win])*100

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B3', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=win, col=0, rotate='top',
							title=freq+' - win'+str(win)+' th'+th)
		#Define a source object with a different color by subject
		xyz = arch_sig['s_xyz']
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_c.visible_obj = False
		sc.add_to_subplot(s_obj_c, row=win, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - BOTTOM VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_bot = BrainObj('B3', translucent=False, hemisphere='both')
		b_obj_bot.alpha = 0.09
		sc.add_to_subplot(b_obj_bot, row=win, col=1, rotate='bottom')
		s_obj_b = SourceObj('bottom', xyz, mask=mask,data=data)
		b_obj_bot.project_sources(s_obj_b,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_b.visible_obj = False
		sc.add_to_subplot(s_obj_b, row=win, col=1)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B3', translucent=False, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=win, col=2, rotate='left',title='Right hemisphere')
		s_obj_r = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r.project_sources(s_obj_r,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=win, col=2)

		b_obj_r2 = BrainObj('B3', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=win, col=3, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_r2.visible_obj = False
		sc.add_to_subplot(s_obj_r2, row=win, col=3)
		
		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B3', translucent=False, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=win, col=4, rotate='left', title='Left hemisphere')
		s_obj_l = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l.project_sources(s_obj_l,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=win, col=4)

		b_obj_l2 = BrainObj('B3', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=win, col=5, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap='jet', clim=(-30,30),radius=12.)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=win, col=5)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='ΔPower (%)', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=win, col=6, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method, min_sig,bsl,freq,th),autocrop=True,print_size=(9,12))