from visbrain import Brain, Colorbarload
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '05', 'Locpow', True
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_E_EpiPerf_LowHigh/')
path_mask = join(path_npz, 'masks_visbrain'+th+'/')
path_to_save = join(path_npz, 'Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}_corr.npy')
f_form_save = join(path_to_save, '{}_sources_min{}_win{}_'+feat+'_allfreqs_th{}_corr.png')
###############################################################################

freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
wins = [3]
min_sigs = [2]
methods = ['s_Mai_RL'] #'aal_RL','BA'

for win,min_sig, method in product(wins,min_sigs, methods):
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
		arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
		elecs = arch_sig['s_elec']
		if exists(masks_vis_form.format('MAI_RL',min_sig,freq,str(win),th)):
			mask = np.load(masks_vis_form.format('MAI_RL',min_sig,freq,str(win),th))
		else:
			mask = np.ones(elecs.shape[0])#create a mask for all elecs when non significant
		#print(mask)

		# =============================================================================
		#						Electrodes sources by subject - RIGHT HEM - EXT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='bottom',
							title=freq[2:]+' win '+str(win)+' th'+th,title_size=10)
		#Define a source object with a different color by subject
		xyz, subjects,channels = arch_sig['s_xyz'],arch_sig['su_codes'], arch_sig['s_channels']
		xyz_sig = xyz[np.where(mask==False)]
		subjects_sig = subjects[np.where(mask==False)]
		elecs_labels_sig=None
		u_color = ["darkblue", "royalblue", "deepskyblue", "mediumspringgreen", "yellow", "darkorange","red"]
		color = [u_color[int(k[1])] for k in subjects_sig]
		data = np.arange(len(subjects_sig))
		s_obj_c = SourceObj('Color', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data, text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		sc.add_to_subplot(s_obj_c, row=i, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - LEFT HEM - EXT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=1, rotate='left')
		s_obj_r = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_r.set_visible_sources('left')
		sc.add_to_subplot(s_obj_r, row=i, col=1)

		# =============================================================================
		#						Electrodes sources by subject - RIGHT HEM - EXT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=i, col=2, rotate='right')
		s_obj_r = SourceObj('S_right', xyz_sig, symbol='disc', color=color,radius_min=10., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
		s_obj_r.set_visible_sources('right')
		sc.add_to_subplot(s_obj_r, row=i, col=2)

	#sc.preview()
sc.screenshot(f_form_save.format(method, min_sig,win,th),autocrop=True,print_size=(4,10))