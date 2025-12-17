from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import BrainObj,SourceObj, ConnectObj
from visbrain.io import download_file, path_to_visbrain_data

feat = 'pow'
###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'Visbrain/npz_files_3win_expi/')
path_to_save = join(path, 'Visbrain/All_subjects/18rois_3wins/')
f_form = '{}_sources_{}_odor_EpiScore_Expi_phys_None.npz'
f_form = join(path_npz, f_form)
###############################################################################

# load all elec coordinates to initiate Brain
hem = 'both'
arch_sig = np.load(f_form.format('All_subjects','0_VLFC'))
xyz = arch_sig['s_xyz']
subjects = arch_sig['su_codes']
b_obj = BrainObj('B2', hemisphere=hem, translucent=True)
b_obj.alpha = 0.09

# Initiate Brain with xyz electrode coordinates
u_color = ["darkblue", "royalblue", "deepskyblue", "mediumspringgreen", "yellow", "darkorange","red"]
color = [u_color[int(k[1])] for k in subjects]
data_plot = np.arange(len(subjects))
s_obj = SourceObj('S_obj', xyz, edge_color='white', symbol='disc', 
					edge_width=.5, radius_min=12., alpha=.9,text_size=1.5, text_color='white')
s_obj.visible_obj = True
#s_obj.text = arch_sig['s_elec_num_su']
vb = Brain(brain_obj=b_obj, source_obj=[s_obj]) #project_mask_color='gray'

# Parameters to update the brain'
freqs = ['0_VLFC', '1_delta', '2_theta', '3_alpha', '4_beta','5_gamma1','6_gamma2']
phase = 'all'

# Update which electrodes to plot and sources data
for freq in freqs:
	# Load the npz file with su_codes, xyz of signif and not signif elec
	mat = np.load(f_form.format('All_subjects',freq))
	#define a DA mask vector to control which sources are displayed and their colors
	s_th_da, da = mat['s_th_05'], mat['s_da']
	mask = np.array([])
	for s in range(7):
		n_elecs = s_th_da[np.where(subjects=='S'+str(s))].shape[0]
		print (s, n_elecs)
		for i in range(n_elecs): # loop across electrodes
			da_sig = np.ravel(da[np.where(subjects=='S'+str(s))][i])
			th_da = max(s_th_da)
			underp = np.where(da_sig > th_da)[0]
			pvsplit = np.split(underp, np.where(np.diff(underp) != 1)[0]+1)
			signif = [False if len(k) >= 2 else True for k in pvsplit ]
			signif = [False if False in signif else True]
			mask = np.hstack((mask, signif)) if np.size(mask) else signif
	vb.sources_control('S_obj', data=data_plot, mask=mask.astype(bool), mask_color='gray', 
		visible=True, color=color)
	df = s_obj.analyse_sources(roi_obj='aal',distance=15., keep_only=['Amygdala (R)', 'Cingulum Ant (R)', 
	'Cingulum Mid (R)', 'Frontal Inf Orb (R)', 'Frontal Inf Tri (R)', 'Frontal Mid (R)', 'Frontal Mid Orb (R)', 
	'Frontal Sup (R)','Frontal Sup Medial (R)', 'Frontal Sup Orb (R)', 'Hippocampus (R)', 'Insula (R)','ParaHippocampal (R)', 
	'Temporal Inf (R)', 'Temporal Mid (R)', 'Temporal Pole Mid (R)','Temporal Pole Sup (R)', 'Temporal Sup (R)',
	'Amygdala (L)', 'Cingulum Ant (L)', 'Cingulum Mid (L)', 'Frontal Inf Orb (L)', 'Frontal Inf Tri (L)', 
	'Frontal Mid (L)', 'Frontal Mid Orb (L)', 'Frontal Sup (L)','Frontal Sup Medial (L)', 'Frontal Sup Orb (L)', 
	'Hippocampus (L)', 'Insula (L)','ParaHippocampal (L)','Temporal Inf (L)', 'Temporal Mid (L)', 
	'Temporal Pole Mid (L)','Temporal Pole Sup (L)', 'Temporal Sup (L)'])
	print(len(df))

	rots =  {'left':'left','right':'right','top':'all'}
	for rot in rots:
		print(rot,rots[rot])
		b_obj.rotate(fixed=rot)
		s_obj.set_visible_sources(rots[rot])
		vb.screenshot(f_form_save.format('All_subjects',freq,phase,rots[rot],rot), transparent=False,
              autocrop=True, print_size=(3, 3), unit='inch',dpi=200.,
              bgcolor='black', canvas='main')

