from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

score, feat = 'Epi', 'AUC'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_odor_expi_AAL_selection/')
path_mask = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score_5aal/masks_visbrain/')
path_to_save = join(path, 'Visbrain/All_subjects/1_New_allrois_time_'+score+'Score_5aal/Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_'+score+'Score_Expi_sel_aal_5_{}.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_'+feat+'_bsl_{}_{}_3wins_th{}.png')
###############################################################################

freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha','4_beta','5_gamma1','6_gamma2'] #
bsls = ['None']
wins = [3,4,5]
ths = ['05']
min_sigs = [0] #for min=0 all methods are equivalent
methods = ['aal_RL']#'aal','aal_RL',BA
clim_min, clim_max = 0.5, 1.
cmap = 'autumn'

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
	for i,win in enumerate(wins):
		mask = np.load(masks_vis_form.format(method,min_sig,freq,str(win),th))
		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B3', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='top',
							title=freq+' - win'+str(win)+' th'+th)
		#Define a source object with a different color by subject
		xyz, da = arch_sig['s_xyz'], arch_sig['s_da']
		data = [max(da[i]) for i in range(da.shape[0])]
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
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
		b_obj_bot.project_sources(s_obj_b,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
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
		b_obj_r.project_sources(s_obj_r,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=2)

		b_obj_r2 = BrainObj('B3', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=3, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
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
		b_obj_l.project_sources(s_obj_l,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=4)

		b_obj_l2 = BrainObj('B3', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=5, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap=cmap, clim=(clim_min, clim_max),radius=12.,
			vmin=0.6,under='grey')
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=i, col=5)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='Max AUC Score', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=6, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method, min_sig,bsl,freq,th),autocrop=True,print_size=(9,12))