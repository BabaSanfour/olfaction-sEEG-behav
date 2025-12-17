"""
Compute distances between trials (space, time, odor, respi)
"""
import numpy as np
from os import path
import pandas as pd
from itertools import combinations

from brainpipe.system import study
from params import subjects
from codes import num_odor_to_ref
from utils import pos_cirles, odor_su_score

def reorder_df(df_to_reorder, od_order):
    """
    order : array of indexes (int)
    """
    #order df according to mat file
    df_to_reorder['order'] = [np.where(od_order==od)[0][0] for od in df_to_reorder['odor_num']]
    df_reordered = df_to_reorder.sort_values(by=['order','trial_index'])
    return df_reordered

def compute_distance2D(df,r1,r2,col0,col1):
    xy1 = np.array([df[col0].iloc[r1].values[0],df[col1].iloc[r1].values[0]])
    xy2 = np.array([df[col0].iloc[r2].values[0],df[col1].iloc[r2].values[0]])
    dist = np.round(np.linalg.norm(xy1-xy2),2)
    return dist

def compute_distance1D(df,r1,r2,col0):
    x1 = np.array([df[col0].iloc[r1].values[0]])
    x2 = np.array([df[col0].iloc[r2].values[0]])
    dist = np.round(np.linalg.norm(x1-x2),2)
    return dist

def compute_dist_trials_su_E(df_perf, df_od_dist, su, cols, norm=True):
    """
    Compute trials distance in terms of odors properties, breathing, space and time

    Parameters
    ----------
    df_perf : data by trials behavior
    df_od_list : dist for all pairs of odors
                sum(intensity, pleasantness, familiarity, type, descriptors)
    su : subject ID
    cols : columns to take to compute distance (in behav df)
    norm : whether or not rescale distance between 0 and 1

    Returns
    -------
    df_save : df with all distance for all targeted trials
    """
    df_E = df_perf
    df_E['0_insp_V'] = np.abs(df_E['total_encodage_0_insp_volume'])
    df_E['0_exp_V'] = df_E['total_encodage_0_exp_volume']*-1
    df_E = df_E[cols]
    df_E = df_E.fillna(df_E.mean())
    df_od = pd.read_csv(df_od_dist)

    pl_d, fam_d, int_d, stru_d, type_d, od_dist = [], [], [], [], [], []
    od_pairs, temp_dist, insp_dist, resp_dist = [], [], [], []
    spa_dist, rich_dist = [], []
    for r1, r2 in combinations(np.arange(df_E.shape[0]),2):
        o1 = str(df_E['odor_num'].iloc[r1])
        o2 = str(df_E['odor_num'].iloc[r2])
        p1 = 'P'+str(num_odor_to_ref[su][int(o1)][1])
        p2 = 'P'+str(num_odor_to_ref[su][int(o2)][1])

        if o1 != o2:
            comb, comb_inv = o1+'_'+o2, o2+'_'+o1
            od_pairs.append(comb)
            pl_d.append(df_od['Pl_d'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])
            fam_d.append(df_od['Fam_d'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])
            int_d.append(df_od['Int_d'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])
            stru_d.append(df_od['Stru_d'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])
            type_d.append(df_od['Type_d'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])
            od_dist.append(df_od['dist_sum'].loc[df_od['od_pairs'].isin([comb,comb_inv])].values[0])

            temp_dist.append(compute_distance2D(df_E,r1,r2,
                                              col0=['trial_time'],col1=['encoding_day']))
            insp_dist.append(compute_distance1D(df_E,r1,r2,col0=['0_insp_V']))
            resp_dist.append(compute_distance2D(df_E,r1,r2,
                                              col0=['0_insp_V'],col1=['0_exp_V']))
            coord1, coord2 = np.array(pos_cirles[p1]), np.array(pos_cirles[p2])
            spa_dist.append(np.round(np.linalg.norm(coord1-coord2),2))
            score1 = np.array(odor_su_score[su][int(o1)])
            score2 = np.array(odor_su_score[su][int(o2)])
            rich_dist.append(np.round(np.linalg.norm(score1-score2),2))

    data = np.concatenate((np.array(od_pairs)[:,np.newaxis],
                           np.array(pl_d)[:,np.newaxis],
                           np.array(fam_d)[:,np.newaxis],
                           np.array(int_d)[:,np.newaxis],
                           np.array(stru_d)[:,np.newaxis],
                           np.array(type_d)[:,np.newaxis],
                           np.array(od_dist)[:,np.newaxis],
                           np.array(temp_dist)[:,np.newaxis],
                           np.array(resp_dist)[:,np.newaxis],
                           np.array(spa_dist)[:,np.newaxis],
                           np.array(rich_dist)[:,np.newaxis]),axis=1)
    df_dist = pd.DataFrame(data, columns=['od_pairs','Pl_d','Fam_d','Int_d',
                              'Stru_d','Type_d','od_dist','temp_dist',
                              'resp_dist','spa_dist','rich_dist'])
    if norm == True:
        #rescale all distances to be between 0 and 1
        for c in df_dist.columns:
            if c != 'od_pairs':
                df_dist[c] = df_dist[c].astype(float)
                df_dist[c] = [(x - min(df_dist[c]))/(max(df_dist[c]-min(df_dist[c]))) \
                                for x in df_dist[c]]

    return df_dist


if __name__ == "__main__":

    subjects = ['CHAF','FERJ','LEFC','PIRJ','SEMC','VACJ']
    PATH = '/media/1_Analyses_Intra_EM_Odor/1bis_OE_BaseSam/JPlailly201306_seeg_ALS/behavior/'
    df_perf = 'encoding_individual_results.xls'
    df_od_dist = 'distance_odors_all_pairs_meth=sum.csv'
    df_save = 'distance_dims=all_odors=all_su={}_norm={}.csv'

    for su in subjects:

        df_all_dist = compute_dist_trials_su_E(df_perf, df_od_dist, su, trials_sel)
