"""Colorer des ROI en fonction de données qui leur sont attribuées."""
import numpy as np
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt
from visbrain.objects import SceneObj, BrainObj, RoiObj, ColorbarObj
from visbrain.utils import array2colormap
from os.path import exists

#=============================================================================
path = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/Bilan_classif/'
path_rois = path + '1_signif_elecs_wins_phase/'
path2save = path+ '3_Figures_Phase/'
#=============================================================================
conds = {'1':['poor','partial',1,2,3,4,5],'2':['partial','detailed',1,2,3,4,5],
         '3':['poor','detailed',1,2,3,4,5]}
freqs = ['2_theta'] #'2_theta','6_gamma2','0_VLFC','1_delta','3_alpha'
rots = ['front','left','right']
#=============================================================================

for freq in freqs:
	# =============================Create a default scene ========================
	CAM_STATE = dict(azimuth=0,        # azimuth angle
	                 elevation=90)     # elevation angle
	CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
	                  rect=(-.3, -2., 1., 4.))
	sc = SceneObj(camera_state=CAM_STATE, bgcolor=(.1, .1, .1),size=(1000, 750)) #largeur, hauteur

	for i,cond in enumerate(conds):
		for j,rot in enumerate(rots):
			# ============================Create a brain object ============================
			b_obj = BrainObj('B1', translucent=True, hemisphere='both')
			b_obj.alpha = 0.09
			sc.add_to_subplot(b_obj, row=j, col=i, rotate=rot)#title="ROIs decoding "+conds[cond][0]+' vs '+conds[cond][1])

			# ============================Create a roi object ==============================
			r_obj = RoiObj('aal')
			labels = r_obj.get_labels() #df of all rois in the atlas

			#load selected rois to plot
			roisname = cond+'_Classif_'+conds[cond][0]+'_'+conds[cond][1]+'_'+freq+'_win1.csv'
			roisfile = pd.read_csv(path_rois+roisname)
			grouped = roisfile.groupby(['s_aal_RL']).nunique()

			#select labels to plot
			lim = conds[cond][2]
			df_sel = grouped.query('su_codes >= @lim')
			idx_to_plot, data_to_plot = df_sel.index.values, df_sel['su_codes'].values
			idx_to_plot = [e for e in idx_to_plot if e != 'Not f']
			data_to_plot = np.repeat(data_to_plot,2)
			idx_to_plot = [s+suf for s,suf in product(idx_to_plot,[' (R)',' (L)'])]
			print('rois to plot', idx_to_plot)
			index = labels[labels['aal'].isin(idx_to_plot)].index.values
			print('id of rois',index)
				
			#convert data to plot in colormap
			roi_color = array2colormap(data_to_plot, cmap='plasma', 
				clim=(conds[cond][2],conds[cond][5]))
			
			#create dict of colors
			dico_color = {k: i for k, i in zip(index, roi_color)}
			r_obj.select_roi(select=index, roi_to_color=dico_color, smooth=15)
			sc.add_to_subplot(r_obj, row=j, col=i)

		# # =============================================================================
		# #						ColorBAR for AUC modulation 
		# # =============================================================================
		#cb_rep = ColorbarObj(r_obj, cblabel='Nb of pqtients sig for Phase', border=False, **CBAR_STATE)
		#sc.add_to_subplot(cb_rep, row=j, col=i+1, width_max=200, height_max=600)

		#sc.preview()
	sc.screenshot(path2save+'AAL_ROIs_Phase_win1_'+freq[2:]+'.png',autocrop=True)
		
