from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th, bsl, feat = '01', '700ms', 'AUC'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_E_EpiPerf_LowHigh_'+bsl+'/')
path_mask = join(path_npz, 'masks_visbrain'+th+'/')
path_to_save = join(path_npz, 'Brain_Plots/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT_{}.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_'+feat+'_bsl_{}_th{}.png')
###############################################################################

freqs = ['2_theta','6_gamma2'] #'0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2'
list_su = ['S0','S1','S2','S3','S4','S5','S6']
wins, radius = [1], 12.
min_sigs = [1]
methods = ['s_Mail_RL']

for freq,min_sig, method in product(freqs,min_sigs, methods):
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
	arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor',bsl))
	#print(arch_sig.files)
	subjects = arch_sig['su_codes']
	#find the lowest significant perm threshold
	perm_all = np.array([])
	for su in list_su:
	    #print(arch_sig.files)
	    s_th = arch_sig['s_th_'+th][np.where(subjects==su)]
	    perm_all = np.hstack((perm_all, max(s_th))) if np.size(perm_all) else max(s_th)
	perm_min = round(min(perm_all),1)

	for win in wins:
		mask = np.load(masks_vis_form.format('MAI_RL',min_sig,freq,str(win),th))
		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=win, col=0, rotate='bottom',
							title=freq+' - win'+str(win)+' th'+th)
		#Define a source object with a different color by subject
		xyz, idx = arch_sig['s_xyz'], np.where(np.max(arch_sig['s_da'], axis=1))#[:,win])
		data = idx#arch_sig['s_da'][:,idx]
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap='autumn', clim=(0.5,1.),
			vmin=perm_min, under='grey',radius=radius)
		s_obj_c.visible_obj = False
		sc.add_to_subplot(s_obj_c, row=win, col=0)
		#sc.preview()

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_r = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r.alpha = 0.09
		sc.add_to_subplot(b_obj_r, row=win, col=1, rotate='left',title='Right hemisphere')
		s_obj_r = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r.project_sources(s_obj_r,'modulation', cmap='autumn', clim=(0.5,1.),
			vmin=perm_min, under='grey',radius=radius)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=win, col=1)

		b_obj_r2 = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=win, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap='autumn', clim=(0.5,1.),
			vmin=perm_min, under='grey',radius=radius)
		s_obj_r2.visible_obj = False
		sc.add_to_subplot(s_obj_r2, row=win, col=2)
		
		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=win, col=3, rotate='left', title='Left hemisphere')
		s_obj_l = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l.project_sources(s_obj_l,'modulation', cmap='autumn', clim=(0.5,1.),
			vmin=perm_min, under='grey',radius=radius)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=win, col=3)

		b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=win, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap='autumn', clim=(0.5,1.),
			vmin=perm_min, under='grey',radius=radius)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=win, col=4)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='AUC Score', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=win, col=5, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format('All_subjects',min_sig,freq,th),autocrop=True,print_size=(9,12))