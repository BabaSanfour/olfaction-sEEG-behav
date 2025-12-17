from visbrain import Brain, Colorbar
from os import listdir
from brainpipe.statistics import *
from os.path import isfile, join, exists
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

th,thresh, feat, corr = '0.01',0.01, 'pow', False
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/0_Classif_Power_R_EpiPerf_LowHigh_1000perm/')
path_mask = join(path_npz, 'masks_stat/')
path_to_save = join(path_npz, 'Brain_Plots'+th+'/')
npz_form = join(path_npz, '{}_sources_{}_{}_low_high_sel_physFT.npz')
masks_form = join(path_mask, 'All_subjects_mask_stat_{}_minwin{}_th{}.npy') 
f_form_save = join(path_to_save, '{}_proj_intime_min{}_win{}_'+feat+'_{}_th{}.png') 
###############################################################################

freqs = ['2_theta','3_alpha','4_beta','5_gamma1','6_gamma2']
wins, radius = [1], 10.
min_sigs, step = [1], 5#step to decide how many time points to plot
methods = ['s_Mai_RL']
clim = (-100,100)
#times = ['250','750','1250','1750','2225','2750']
times = ['0','500','1000','1500','2000']

#Determine the max permutation threshold
perm_max = np.array([])
for freq in freqs:
    mat = np.load(npz_form.format('All_subjects',freq,'odor'))
    s_perm = mat['s_perm']#[:,:,5:30]
    for elec in range(s_perm.shape[0]):
        perm = s_perm[elec]
        th_perm = perm_pvalue2level(perm, p=thresh, maxst=True)
        perm_max = np.hstack((perm_max,th_perm[0])) if np.size(perm_max) else th_perm[0]
perm_thr = max(perm_max) #find the max of perm thr across time and subjects
print('perm max size',perm_max.shape,'thr max',perm_thr)

for win, freq, min_sig, method in product(wins,freqs,min_sigs,methods):
	arch_sig = np.load(npz_form.format('All_subjects',freq, 'odor'))
	#print(arch_sig.files)
	subjects = arch_sig['su_codes']
	power = ((arch_sig['s_elec_pow1'][:,:]-arch_sig['s_elec_pow0'][:,:])/abs(arch_sig['s_elec_pow0'][:,:]))*100
	da = arch_sig['s_da'][:,:]#[5:30]
	print('da shape',da.shape)
	
	# =============================================================================
	#								Create a default scene
	# =============================================================================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90,     # elevation angle
	                 )
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(700, 800)) #largeur, hauteur
	# =============================================================================

	for i,t in enumerate(np.arange(0,da.shape[1],step)): #loop across time for each time point
		print(i,t)
		mask,da_plot,pow_plot = [], np.array([]),np.array([])
		for elec in range(da.shape[0]):
			da_win = da[elec,t:t+step]
			pow_win = power[elec,t:t+step]
			da_t = max(da_win)
			pow_t = np.mean(pow_win[np.where(da_win==da_t)], axis=0)
			da_plot = np.vstack((da_plot,da_t)) if np.size(da_plot) else da_t
			pow_plot = np.vstack((pow_plot,pow_t)) if np.size(pow_plot) else pow_t
			if da_t > perm_thr:
				mask.append(False)
			else :
				mask.append(True)
		data = pow_plot
		# =============================================================================
		#						Electrodes sources by subject - TOP VIEW 
		# =============================================================================
		#Define a brain object to plot
		b_obj_top = BrainObj('B2', translucent=False, hemisphere='both')
		b_obj_top.alpha = 0.09
		sc.add_to_subplot(b_obj_top, row=i, col=0, rotate='bottom',
							title=freq+' - win'+str(win)+' th'+th+' // '+times[i]+'-'+times[i+1]+'ms')
		#Define a source object with a different color by subject
		xyz = arch_sig['s_xyz']
		s_obj_c = SourceObj('modulation', xyz, data=data,mask=mask)
		b_obj_top.project_sources(s_obj_c,'modulation', cmap='jet', clim=clim,
			radius=radius)
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
		b_obj_r.project_sources(s_obj_r,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_r.visible_obj = False
		sc.add_to_subplot(s_obj_r, row=i, col=1)

		b_obj_r2 = BrainObj('B2', translucent=False, hemisphere='right')
		b_obj_r2.alpha = 0.09
		sc.add_to_subplot(b_obj_r2, row=i, col=2, rotate='right')
		s_obj_r2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_r2.project_sources(s_obj_r2,'modulation', cmap='jet', clim=clim,
			radius=radius)
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
		b_obj_l.project_sources(s_obj_l,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_l.visible_obj = False
		sc.add_to_subplot(s_obj_l, row=i, col=3)

		b_obj_l2 = BrainObj('B2', translucent=False, hemisphere='left')
		b_obj_l2.alpha = 0.09
		sc.add_to_subplot(b_obj_l2, row=i, col=4, rotate='right')
		s_obj_l2 = SourceObj('S_right', xyz, data=data, mask=mask)
		b_obj_l2.project_sources(s_obj_l2,'modulation', cmap='jet', clim=clim,
			radius=radius)
		s_obj_l2.visible_obj = False
		sc.add_to_subplot(s_obj_l2, row=i, col=4)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		cb_rep = ColorbarObj(b_obj_l2, cblabel='Time', **CBAR_STATE)
		sc.add_to_subplot(cb_rep, row=i, col=5, width_max=200, height_max=600)

	#sc.preview()
	sc.screenshot(f_form_save.format(method,min_sig,win,freq,th),autocrop=True,print_size=(9,12))