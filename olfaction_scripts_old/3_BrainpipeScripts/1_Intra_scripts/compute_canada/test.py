# coding: utf-8

import numpy as np
from itertools import product


#PATH = "/home/alsaive/projects/def-kjerbi/alsaive/Intra_EM/"
PATH = "/home/karim/mp2/Intra_EM/"

conds,subjects = ['low','high'],['FERJ','MICP','VACJ','SEMC','LEFC','PIRJ','CHAF']

for su in subjects:
    pow_list = []
    #=========================== Load Power files (nfreq, nelec, nwin, ntrial) =================================    
    mat0 = np.load(PATH+'dataE/'+su+'_odor_'+conds[0]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat0['xpow'][:,:,17:52,:])
    nelecs = mat0['xpow'].shape[1]
    mat1 = np.load(PATH+'dataE/'+su+'_odor_'+conds[1]+'_bipo_sel_physFT_pow.npz')
    pow_list.append(mat1['xpow'][:,:,17:52,:])

    # =========================== Select Power for 1 elec 1 freq =================================                 
    for elec_num, ti in product(range(nelecs),range(35)):
        pow_data_elec = []
        for power in pow_list:
            dat = power[2:7,elec_num,ti,:].swapaxes(0,1)
            pow_data_elec.append(dat)
        x = np.concatenate(pow_data_elec, axis=0)
        y = np.hstack([np.array([i]*power.shape[0]) for i, power in enumerate(pow_data_elec)])
        print(x.shape,y)
        