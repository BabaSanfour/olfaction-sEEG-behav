from visbrain import Brain, Colorbar
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product

from visbrain.objects import SceneObj, BrainObj,SourceObj, ConnectObj, ColorbarObj
from visbrain.io import download_file, path_to_visbrain_data

###############################################################################
path = r'/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Olfacto/'
path_npz = join(path, 'database/Encoding_EpiPerf_4500_expi_noart/')
f_form = '{}_odor_poor_bipo_sel.npz'
f_form = join(path_npz, f_form)
f_form_save = join(path_npz, 'Elecs_{}_labels_{}.csv')
###############################################################################

subjects = ['LEFC','CHAF','VACJ','SEMC','FERJ','MICP','PIRJ']

for su in subjects:
	arch_sig = np.load(f_form.format(su))
	print(arch_sig.files)
	xyz, labels, channels = arch_sig['xyz'], arch_sig['label'], arch_sig['channel']
	s_obj_c = SourceObj('Color', xyz)
	meth = 'brodmann'
	df = s_obj_c.analyse_sources(roi_obj=meth,distance=15.)
	df['labels'], df['channels'] = labels, channels
	df.to_csv(f_form_save.format(su,meth))