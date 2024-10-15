from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product
import pandas as pd

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data


##############################################################################
roi, freq, exp, pval = 'OFC_olf', 'theta', 'E', '0.01'
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'figure/TPSim_LDA_'+exp+'_by_cond_6freqs_3s_dissim/npy_figs/')
path_to_save = join(path_npz, 'da_plots/')
f_form = join(path_npz, '{}_{}_{}_{}_tps.csv'.format(freq,roi,exp,pval))
f_form_save = join(path_to_save, roi+'_'+freq+'_'+exp+'_plot_da_sources.png')
cmap = 'autumn'
clim = (0.6,0.8)
##############################################################################

# =============================================================================
#								Create a default scene
# =============================================================================
CAM_STATE = dict(azimuth=0,        # azimuth angle
                 elevation=90,     # elevation angle
                 )
CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                  rect=(-.3, -2., 1., 4.))
sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1,1,1),size=(1100, 400))

df = pd.read_csv(f_form)
x,y,z = df['x'].values[:,np.newaxis], df['y'].values[:,np.newaxis], df['z'].values[:,np.newaxis]
xyz = np.concatenate((x,y,z),axis=1)
data = df[freq+'_AUC'].values
print(data)

# =============================================================================
#						Electrodes sources by subject - TOP VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_top = BrainObj('B2', translucent=True)
b_obj_top.alpha = 0.09
sc.add_to_subplot(b_obj_top, row=0, col=0, 
					title='Electrodes localization by subject - Top')
#Define a source object with a different color by subject
labels = None
s_obj_c = SourceObj('modulation', xyz, symbol='disc',radius_min=10., radius_max=11., 
			alpha=.95, text=labels, text_size=12, text_color='white' )
s_obj_c.color_sources(data=data, cmap=cmap, clim=clim)
sc.add_to_subplot(s_obj_c, row=0, col=0)
#sc.preview()

# =============================================================================
#						Electrodes sources by subject - RIGHT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
b_obj_r.alpha = 0.09
sc.add_to_subplot(b_obj_r, row=0, col=1, rotate='right',
					title='Electrodes localization by subject - Right')
s_obj_r = SourceObj('S_right', xyz, symbol='disc', radius_min=10., 
			alpha=.95, data=data, text=labels)
s_obj_r.set_visible_sources('right')
s_obj_r.color_sources(data=data, cmap=cmap, clim=clim)
sc.add_to_subplot(s_obj_r, row=0, col=1)

# =============================================================================
#						Electrodes sources by subject - LEFT VIEW 
# =============================================================================
#Define a brain object to plot
b_obj_l = BrainObj('B2', translucent=True, hemisphere='left')
b_obj_l.alpha = 0.09
sc.add_to_subplot(b_obj_l, row=0, col=2, rotate='left',
					title='Electrodes localization by subject - Left')
s_obj_l = SourceObj('S_right', xyz, symbol='disc', radius_min=10., 
			alpha=.95, data=data, text=labels)
s_obj_l.set_visible_sources('left')
s_obj_l.color_sources(data=data, cmap=cmap, clim=clim)
sc.add_to_subplot(s_obj_l, row=0, col=2)

#sc.preview()
sc.screenshot(f_form_save,autocrop=True,print_size=(6,7))
