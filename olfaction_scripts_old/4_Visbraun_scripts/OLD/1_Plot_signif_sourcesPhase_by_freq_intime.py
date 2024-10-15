from visbrain import Brain, Colorbar
from os import listdir, makedirs
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

conds = ['good', 'bad']
feat = 'PHASE'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figures_npz/1_'+conds[0]+'_'+conds[1]+'_4500_expi_noart/')
path_mask = join(path_npz, 'masks_visbrain/')
path_to_save = join(path_npz, 'Visbrain_Plots/')
npz_form = join(path_npz,'{}_sources_phase_{}_odor_'+conds[0]+'_'+conds[1]+'_Expi_sel_phys.npz')
masks_vis_form = join(path_mask, '{}_mask_min{}_{}_minwin{}_th{}.npy')
f_form_save = join(path_to_save, '{}_signif_min{}_'+feat+'_bsl_{}_{}_th{}_'+conds[0]+'_'+conds[1]+'.png')
###############################################################################
if not exists(path_to_save):
	makedirs(path_to_save)
###############################################################################

freqs = ['0_VLFC','1_delta','2_theta', '3_alpha']
bsls = ['None']
wins = [1,2,3]
ths = ['05']
min_sigs = [2,3,4,5]
methods = ['aal_RL','BA']#'aal','aal_RL',BA
clim_min, clim_max = 0.5,1
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
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(1000, 500)) #largeur, hauteur
	# ============================================================================
	arch_sig = np.load(npz_form.format('All_subjects',freq))
	#print(arch_sig.files)

	for i,win in enumerate(wins):
		#load data to plot and mask
		mask = np.load(masks_vis_form.format(method,min_sig,freq,str(win),th))
		print(mask.dtype, len(mask))
		if mask.dtype == bool:
			da = arch_sig['s_da']
			da_max, idx_all = np.array([]), np.array([])
			for elec in range(da.shape[0]):
				da_elec = max(da[elec])
				idx = [i for i,j in enumerate(da[elec]) if j ==da_elec]
				da_max = np.hstack((da_max, da_elec)) if np.size(da_max) else da_elec
				idx_all = np.hstack((idx_all, np.mean(idx))) if np.size(idx_all) else idx
			data = da_max[np.where(mask==False)]
			
			# =============================================================================
			#						Electrodes sources by subject - TOP VIEW 
			# =============================================================================
			#Define a brain object to plot
			b_obj_top = BrainObj('B3', translucent=True, hemisphere='both')
			b_obj_top.alpha = 0.09
			sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='top',
								title=freq+' - minwin'+str(win)+' th'+th)
			# #Define a source object with a different color by subject
			xyz = arch_sig['s_xyz']
			xyz_sig = xyz[np.where(mask==False)]
			elecs_labels_sig=None
			s_obj_c = SourceObj('modulation', xyz_sig, radius_min=10.,radius_max=11.)
			s_obj_c.color_sources(data=data, cmap=cmap, clim=(clim_min,clim_max))
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
			s_obj_r = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10.,radius_max=11.,
						alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
			s_obj_r.set_visible_sources('right')
			s_obj_r.color_sources(data=data, cmap=cmap, clim=(clim_min,clim_max))
			sc.add_to_subplot(s_obj_r, row=i, col=1)
			#sc.preview()

			b_obj_r2 = BrainObj('B2', translucent=True, hemisphere='right')
			b_obj_r2.alpha = 0.09
			sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
			s_obj_r2 = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10.,radius_max=11.,
						alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
			s_obj_r2.set_visible_sources('right')
			s_obj_r2.color_sources(data=data, cmap=cmap, clim=(clim_min,clim_max))
			sc.add_to_subplot(s_obj_r2, row=i, col=2)

			# # =============================================================================
			# #						Electrodes sources by subject - LEFT VIEW 
			# # =============================================================================
			#Define a brain object to plot
			b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
			b_obj_l.alpha = 0.09
			sc.add_to_subplot(b_obj_l, row=i, col=3, rotate='left',
								title=freq+' win '+str(win)+' th'+th+' - Left')
			s_obj_l = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10.,radius_max=11., 
						alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
			s_obj_l.set_visible_sources('left')
			s_obj_l.color_sources(data=data, cmap=cmap, clim=(clim_min,clim_max))
			sc.add_to_subplot(s_obj_l, row=i, col=3)

			b_obj_l2 = BrainObj('B2', translucent=True, hemisphere='left')
			b_obj_l2.alpha = 0.09
			sc.add_to_subplot(b_obj_l2, row=i, col=4, rotate='right')
			s_obj_l2 = SourceObj('S_right', xyz_sig, symbol='disc', radius_min=10., radius_max=11.,
						alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,text_color='white', text_size=9)
			s_obj_l2.set_visible_sources('left')
			s_obj_l2.color_sources(data=data, cmap=cmap, clim=(clim_min,clim_max))
			sc.add_to_subplot(s_obj_l2, row=i, col=4)

			# # =============================================================================
			# #						ColorBAR for AUC modulation 
			# # =============================================================================
			cb_rep = ColorbarObj(s_obj_l2, cblabel=feat+'(AUC Score)', border=False, **CBAR_STATE)
			sc.add_to_subplot(cb_rep, row=i, col=5, width_max=200, height_max=600)
		else:
			continue
		#sc.preview()
		sc.screenshot(f_form_save.format(method, min_sig,bsl,freq,th),autocrop=True,print_size=(9,12))
