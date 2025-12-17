from visbrain.gui import Brain
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'database/Encoding_By_Odor/try_PIRJ/')
f_form = '{}_odor_{}_bipo_all_noWM.npz'
f_form = join(path_npz, f_form)
f_form_save = join(path, 'database/Encoding_By_Odor/_labels/Elecs_{}_labels_{}.csv')
###############################################################################

#Download the MIST atlas
#download_file('MIST_ROI.zip', unzip=True, astype='example_data')

#subjects = ['CHAF','LEFC','FERJ','PIRJ','SEMC','VACJ','MICP'] #MICP
subjects = ['PIRJ'] #MICP
#conds = ['high'] #same electrodes in E1 and E2
odors_su = {'CHAF': (5,7,8,9,1,2,3,4),
            'LEFC':(1,2,3,4,14,15,16,17),
            'PIRJ':(4,9,1,18,6,5,7), #missing odor 15
            'VACJ':(14,15,16,17,10,11,12,13),
            'SEMC':(10,11,12,13,5,7,8,9),
            'MICP':(2,12,6,8,3,18,9,14),
            'FERJ':(16,17,5,7,12,13,2,1)}

methods = ['brodmann','aal','talairach']

for su, meth in product(subjects, methods):
	for od in odors_su[su]:
		print (su, od)
		arch_sig = np.load(f_form.format(su,od))
		print(arch_sig.files)
		xyz, labels, channels = arch_sig['xyz'], arch_sig['label'], arch_sig['channel']
		s_obj_c = SourceObj('Color', xyz)
		df = s_obj_c.analyse_sources(roi_obj=meth,distance=15.)
		df['labels'], df['channels'] = labels, channels
		df.to_csv(f_form_save.format(su,meth))
