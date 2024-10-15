from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_3win_expi/')
path_to_save = join(path, 'Visbrain/All_subjects/18rois_3wins/')
f_form = '{}_sources_{}_{}_EpiScore_Expi_phys_{}.npz'
f_form = join(path_npz, f_form)
f_form_save = join(path_to_save, '{}_signif_elecs_bsl_{}_{}_3wins_th{}.png')
mask_form = '{}/{}_mask_stat_{}_win{}_th{}.npy'
mask_form = join(path_to_save, mask_form)
###############################################################################

freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
bsls = ['None','mean_trial']
wins = [0,1,2]
ths = ['05']

for freq, bsl, th in product(freqs,bsls,ths):
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(1000, 700)) #largeur, hauteur

	for win in wins:
		mask = np.load(mask_form.format(bsl,'All_subjects',freq,str(win),th))

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=win, col=0, rotate='top',
							title=freq+' win '+str(win)+' th'+th+' - top')
		#Define a source object with a different color by subject
		arch_sig = np.load(f_form.format('All_subjects',freq, 'odor','None'))
		xyz, subjects = arch_sig['s_xyz'],arch_sig['su_codes']
		u_color = ["darkblue", "royalblue", "deepskyblue", "mediumspringgreen", "yellow", "darkorange","red"]
		color = [u_color[int(k[1])] for k in subjects]
		data = np.arange(len(subjects))
		s_obj_c = SourceObj('Color', xyz, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, mask=mask)
		sc.add_to_subplot(s_obj_c, row=win, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=win, col=1, rotate='left',
							title=freq+' win '+str(win)+' th'+th+' - Right')
		s_obj_r = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, mask=mask)
		s_obj_r.set_visible_sources('right')
		sc.add_to_subplot(s_obj_r, row=win, col=1)

		b_obj_r2 = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=win, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, mask=mask)
		s_obj_r2.set_visible_sources('right')
		sc.add_to_subplot(s_obj_r2, row=win, col=2)

		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=win, col=3, rotate='left',
							title=freq+' win '+str(win)+' th'+th+' - Left')
		s_obj_l = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, mask=mask)
		s_obj_l.set_visible_sources('left')
		sc.add_to_subplot(s_obj_l, row=win, col=3)

		b_obj_l2 = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=win, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, mask=mask)
		s_obj_l2.set_visible_sources('left')
		sc.add_to_subplot(s_obj_l2, row=win, col=4)

	#sc.preview()
	sc.screenshot(f_form_save.format('All_subjects',bsl,freq,th),autocrop=True,print_size=(9,12))