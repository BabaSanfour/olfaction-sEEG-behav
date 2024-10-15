from visbrain.gui import Brain
from brainpipe.statistics import *
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th, tresh, feat, corr = '0.01',0.01, 'tpsAUC', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/LDA_TPSim_E_R_all_BBG/')
path_mask = join(path_npz, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
f_form_save = join(path_to_save, '{}_proj_min{}_win{}_'+feat+'_{}_{}_th{}_BGA.png') 
###############################################################################

freqs = ['3_gamma'] #'0_theta','1_alpha','2_beta','3_gamma'
#rois_to_keep = ['ACC','Amg','Amg-PirT','HC','IFG','Ins','MFG','OFC','PHG',
#                'SFG','pPirT']
rois_to_keep = ['IFG']
wins, radius = [1.0], 10.
min_sigs = [1]
methods = ['s_Mai_RL']
clim = (0.5,0.9)
#thresh = np.min([0.82, 0.77, 0.67, 0.60, 0.80, 0.68, 0.82]) #Encoding
#thresh = np.min([[0.82, 0.80, 0.82, 0.76, 0.80, 0.80]]) #Retrieval
#thresh = np.min([0.80,0.80,0.80,0.76,0.78,0.78]) #Resting retrieval
thresh = np.min([0.61, 0.65, 0.56, 0.64, 0.52, 0.68]) #TPSim

for win,min_sig, method in product(wins,min_sigs, methods):
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

	for i,freq in zip(range(3,4),freqs):
	#for i,freq in enumerate(freqs):
		arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
		idx_rois = np.where([roi in rois_to_keep for roi in arch_sig['s_labels']])
		#id_subj = np.where([su in ['S0','S1','S2','S4','S5','S6'] for su in arch_sig['su_codes'][idx_rois]])
		da = arch_sig['s_da'][idx_rois]#[id_subj]
		#da = arch_sig['s_da'][idx_rois,10:30][0]
		print('da shape', da.shape)

		subjects = arch_sig['su_codes'][idx_rois]#[id_subj]
		xyz = arch_sig['s_xyz'][idx_rois]#[id_subj]

		mask = np.load(masks_form.format(freq,str(win),th))[idx_rois]#[id_subj]
		
		#Find the max AUC score for all electrodes WHEN MULTIPLE TIME POINTS
		da_max = np.array([])
		for elec in range(da.shape[0]):
			da_elec = np.max(da[elec])
			da_max = np.hstack((da_max,da_elec)) if np.size(da_max) else da_elec
		data = da_max
		print(data.shape, mask.shape)
		#data = da

		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='bottom',
							title=freq+' - min_patient'+str(min_sig)+' th'+th, title_color='black')
		#Define a source object with a different color by subject
		
		print('shapes', xyz.shape, data.shape, mask.shape)
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap='autumn', clim=clim,
			vmin=thresh, under='grey',radius=radius)
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
		b_obj_r.project_sources(s_obj_r,'modulation', cmap='autumn', clim=clim,
			vmin=thresh, under='grey',radius=radius)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=1)

		b_obj_r2 = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap='autumn', clim=clim,
			vmin=thresh, under='grey',radius=radius)
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
		b_obj_l.project_sources(s_obj_l,'modulation', cmap='autumn', clim=clim,
			vmin=thresh, under='grey',radius=radius)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=3)

		b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap='autumn', clim=clim,
			vmin=thresh, under='grey',radius=radius)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=i, col=4)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='AUC Score',txtcolor='black', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=5, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method,min_sig,win,freqs[0],rois_to_keep[0],th),autocrop=True,print_size=(9,12))