from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj, RoiObj
from visbrain.io import download_file, path_to_visbrain_data

th, feat, corr = '0.01', 'tps', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/LDA_TPSim_E_R_all_BBG_within/')
path_npz2 = join(path, 'figure/LDA_TPSim_E_R_all_BBG/')
path_sel = join(path_npz2, 'Brain_Plots'+th+'/')
path_mask = join(path_npz2, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
#f_form_save = join(path_to_save, '{}_proj_min{}_win{}_'+feat+'_{}_th{}_{}_{}.png') 
f_form_save = join(path_to_save, 'Brain_plots_by_odor_{}_{}_rad{}mm_sel_btw.png') 
###############################################################################

freqs = ['0_theta','1_alpha','2_beta','3_gamma'] if feat == 'pow' else ['1_alpha','2_beta']
rois_to_keep = ['OFC']

#coordinates to select (+/-5mm)
x0,y0,z0 = 20, 30,-16
x1,y1,z1 = -20, 32,-16
rad = 12 #radius of the sphere around coordinates

win, radius = 1.0, 1.
min_sig = 1
method = 's_Mai_RL'


# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(300,300)) #largeur, hauteur
# ============================================================================

meth = 'in' #'out','in' (inside or outside the olfactory OFC)
count = -1
for freq in freqs:
	clim = (-80,80)
	cmap = 'jet'

	arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
	print(arch_sig.files)#,arch_sig['s_labels'][:10],arch_sig['s_MAI_RL'][:10])
	print(arch_sig['s_xyz'].shape)
	idx_rois = np.where([roi in rois_to_keep for roi in arch_sig['s_labels']])[0]
	#labels = arch_sig['s_labels']
	#labels_updated = ['pPirT' if lab in ['Amg-PirT','Amg'] else lab for lab in labels]
	
	#select significant electrodes
	mask = np.load(masks_form.format(freq,str(win),th))[idx_rois]
	mask = ~mask #inverse the boolean to indicate who to plot
	print(mask.shape)

	#select electrodes in the ROIs
	xyz = arch_sig['s_xyz'][idx_rois][mask]
	subjects = arch_sig['su_codes'][idx_rois][mask]
	x,y,z = xyz[:,0], xyz[:,1], xyz[:,2]
	print(freq,'xyz shape',xyz.shape,'x',x.shape)
	sel = np.load(path_sel+'sel_elec_{}_{}_{}.npy'.format(freq,meth,rad))
	# sel = np.zeros((x.shape), dtype=bool)
	# for i in range(x.shape[0]):
	# 	if all([x0-rad<=x[i]<=x0+rad, y0-rad<=y[i]<=y0+rad, z0-rad<=z[i]<=z0+rad]):
	# 		sel[i] = sel[i] if meth == 'out' else 1
	# 		print('IN',subjects[i],'elec sel',x[i],y[i],z[i]) if meth == 'in' else False
	# 	elif all([x1-rad<=x[i]<=x1+rad, y1-rad<=y[i]<=y1+rad, z1-rad<=z[i]<=z1+rad]):
	# 		sel[i] = sel[i] if meth == 'out' else 1
	# 		print('IN',subjects[i],'elec sel',x[i],y[i],z[i]) if meth == 'in' else False
	# 	else:
	# 		sel[i] = 1 if meth == 'out' else sel[i] #not in the ROI, to be excluded from the plot
	# 		print('OUT',subjects[i],'elec sel',x[i],y[i],z[i]) if meth == 'out' else False
	
	# np.save(path_to_save+'sel_elec_{}_{}_{}.npy'.format(freq,meth,rad),sel)
	xyz = xyz[sel]
	pow0 = arch_sig['s_elec_pow0'][idx_rois][mask][sel]
	pow1 = arch_sig['s_elec_pow1'][idx_rois][mask][sel]
	if len(pow0.shape)>1:
		pow0 = np.mean(pow0,axis=1)
		pow1 = np.mean(pow1,axis=1)
	print('pow shape',pow0.shape,pow1.shape)

	data = ((pow1-pow0)/abs(pow0))*100
	#mask = np.asarray([x if power[i] > 0  else True for i,x in enumerate(mask)])
	print(freq, data.shape,xyz.shape,mask.shape)#idx_all.shape)

	if np.size(data):
		count += 1
		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=True, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=count, col=0, rotate='bottom',
					title=freq+' th'+th,title_color='black', title_size=12)
		elecs_labels_sig=None
		s_obj_c = SourceObj('modulation', xyz, radius_min=10.,radius_max=11.,
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,
					text_color='black', text_size=9)
		s_obj_c.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_c, row=count, col=0)

		# # =============================================================================
		# #						Electrodes sources by subject - LEFT VIEW 
		# # =============================================================================
		#Define a brain object to plot
		b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
		b_obj_l.alpha = 0.09
		sc.add_to_subplot(b_obj_l, row=count, col=1, rotate='left')
		s_obj_l = SourceObj('S_right', xyz, symbol='disc', radius_min=10.,radius_max=11., 
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,
					text_color='black', text_size=9)
		s_obj_l.set_visible_sources('left')
		s_obj_l.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_l, row=count, col=1)

		# =============================================================================
		#						Electrodes sources by subject - RIGHT VIEW 
		# =============================================================================
		b_obj_r2 = BrainObj('B2', translucent=True, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=count, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, symbol='disc', radius_min=10.,radius_max=11.,
					alpha=.95, data=data,text=elecs_labels_sig,text_bold=True,
					text_color='black', text_size=9)
		s_obj_r2.set_visible_sources('right')
		s_obj_r2.color_sources(data=data, cmap=cmap, clim=clim)
		sc.add_to_subplot(s_obj_r2, row=count, col=2)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(s_obj_r2, cblabel=feat, border=False, **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=count, col=3, width_max=200, height_max=600)
	#sc.preview()
sc.screenshot(f_form_save.format(rois_to_keep[0],meth,rad),autocrop=True,print_size=(4,10))
