"""
Script to concatenate ALL imported data together
KEEP info for odors order and related memory conditions (High Mid low)
"""
import numpy as np
from os import path
import pandas as pd
from brainpipe.system import study

from params import su_list_od, odor_groups_3wgth, odor_groups_wgth, subjects
from utils_functions import revert_dico

def concat_odors_all_info(path_data, su):
    """
    Concatenate all data imported and bipolarized

    Parameters
    ----------
    su : npz files to select clean channels info


    Returns
    -------
    data to save :  an npz file with all data concatenated of shape (nelecs x n_pts x ntrials).
                    electrodes information ('label', 'channel', 'xyz','new_label')
                    including order of odors & associated mem perf
    """

    x_tot, od_list, od_rich3, od_rich2 = np.array([]), [], [], []
    for od in su_list_od[su]:
        mat = np.load(path_data+f_concat.format(su,od),allow_pickle=True)
        ntrials = mat['x'].shape[-1]
        od_list.extend([od]*ntrials)
        od_rich2.extend([dico_2gr[su][od]]*ntrials)
        od_rich3.extend([dico_3gr[su][od]]*ntrials)
        x_tot = np.concatenate((x_tot,mat['x']),axis=-1) if np.size(x_tot) else mat['x']

    dict_data = {}
    for file in mat.files:
        dict_data[file] = mat[file]
    dict_data['x'], dict_data['odor'] = x_tot, np.asarray(od_list)
    dict_data['EM_2gr'], dict_data['EM_3gr'] = np.asarray(od_rich2), np.asarray(od_rich3)
    #print(x_tot.shape,len(od_list),len(od_rich3),len(od_rich2))

    return dict_data

if __name__ == "__main__":

    dico_2gr = revert_dico(odor_groups_wgth)
    dico_3gr = revert_dico(odor_groups_3wgth)

    reps =  ['R_odors/','R_rec/', 'E_odors/']
    f_concat = '{}_odor_{}_bipo.npz'
    f_save = '{}_cond=ALL_odors=ALL_bipo_INFO.npz'
    st = study('Ripples')

    for rep in reps:
        path_data = path.join(st.path, 'database/'+rep)
        for su in dico_2gr:
            print('--> processing',rep,su)
            dic_save = concat_odors_all_info(path_data, su)
            np.savez(path.join(path_data, f_save.format(su)), **dic_save)
