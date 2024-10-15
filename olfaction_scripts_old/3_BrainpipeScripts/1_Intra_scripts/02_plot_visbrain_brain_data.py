
from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
path = r'/media/1_Analyses_Intra_EM_Odor/Olfacto/'
path_csv = join(path, 'database/Encoding_all_Odors/fg_models_theta/')
csv_form = join(path_csv, 'csv_su=all_roi=OFC_all_f=[2-9].csv')
f_form_save = join(path_csv, 'OFC_plot_feat={}_su=all_x=abs.png')
##############################################################################
feats = ['CF','BW']
cmap = 'cool'
rad = 10
rotation = ['right','bottom']

for feat in feats:
    clim = (0,4) if feat == 'BW' else (3,7)
    # =============================================================================
    #								Create a default scene
    # =============================================================================
    CAM_STATE = dict(azimuth=0,        # azimuth angle
                     elevation=90,     # elevation angle
                     )
    CBAR_STATE = dict(cbtxtsz=12, txtsz=10., width=.1, cbtxtsh=3.,
                      rect=(-.3, -2., 1., 4.))
    sc = SceneObj(camera_state=CAM_STATE, bgcolor=(1., 1., 1.),size=(800, 400))

    # =============================================================================
    #						Electrodes sources by subject - TOP VIEW
    # =============================================================================
    # #Define a source object with a different color by subject
    df = pd.read_csv(csv_form)
    df = df.loc[df['theta_detect']==1]
    x = np.abs(df['x'])
    x,y,z = x[:,np.newaxis], df['y'].values[:,np.newaxis], df['z'].values[:,np.newaxis]
    xyz = np.concatenate((x,y,z),axis=1)
    data = df[feat].values
    print(max(data),min(data),np.mean(data),np.median(data))

    # # =============================================================================
    # #						Electrodes sources by subject - BOTTOM VIEWS
    # # =============================================================================

    #Define a brain object to plot
    b_obj_r = BrainObj('B2', translucent=True, hemisphere='right')
    b_obj_r.alpha = 0.09
    s_obj_r = SourceObj('Rep_left', xyz, symbol='disc',radius_min=10.)
    s_obj_r.color_sources(data=data, cmap=cmap)
    sc.add_to_subplot (b_obj_r, row=0, col=0, rotate=rotation[0])
    sc.add_to_subplot(s_obj_r, row=0, col=0)

    #Define a brain object to plot
    b_obj_b = BrainObj('B2', translucent=True, hemisphere='right')
    b_obj_b.alpha = 0.09
    s_obj_b = SourceObj('Rep_left', xyz, symbol='disc',radius_min=10.,radius_max=10.)
    s_obj_b.color_sources(data=data, cmap=cmap)
    sc.add_to_subplot (b_obj_b, row=0, col=1, rotate=rotation[1])
    sc.add_to_subplot(s_obj_b, row=0, col=1)

    sc.screenshot(f_form_save.format(feat),autocrop=True,print_size=(6,7))
