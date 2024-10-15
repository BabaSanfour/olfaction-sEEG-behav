import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import os

path_study = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/respi/'

subjects = ['LEFC','CHAF','VACJ','SEMC','FERJ','MICP','PIRJ']
sessions = ['E1','E2']
conds = ['odorall','no_odor']

def length_all_subjects():
	all_length = []
	for su, sess in product(subjects,sessions):
	    #Load cycles times
	    df_cycles = pd.read_excel(path_study+'cycles_by_cond/All_cycles_features_E.xls',sheetname=su+'_'+sess)
	    #length = (df_cycles['exp_time']-df_cycles['insp_time']).values
	    length = (df_cycles['insp_duration']+df_cycles['exp_duration']).values
	    all_length.append(length)
	all_length = np.concatenate(all_length, axis =0)
	print(all_length.shape,all_length.min(),all_length.max())

	fig, ax = plt.subplots()
	y,x  = np.histogram(all_length, bins = np.arange(np.min(all_length),np.mean(all_length)*2,0.1))
	print(y.shape,x.shape)
	ax.plot(x[:-1], y)
	ax.axvline(np.mean(all_length), ls = '-', label='mean = {:5.2f}'.format(np.mean((all_length))))
	ax.axvline(np.median(all_length), ls = '--', label='median = {:5.2f}'.format(np.median(all_length)))
	plt.legend(loc=0)
	plt.show()
	plt.clf()
	plt.close()

def length_by_subj():
	fig, ax = plt.subplots()
	for su in subjects:
		all_length = []
		for sess in sessions:
			#Load cycles times
			df_cycles = pd.read_excel(path_study+'cycles_by_cond/All_cycles_features_E.xls',sheetname=su+'_'+sess)
			length = (df_cycles['exp_duration']+df_cycles['insp_duration']).values
			#length = (df_cycles['exp_time']-df_cycles['insp_time']).values
			all_length.append(length)
		all_length = np.concatenate(all_length, axis =0)
		print(all_length.min(),all_length.max())
		y,x  = np.histogram(all_length, bins = np.arange(np.min(all_length),np.mean(all_length)*2,0.3))
		print(y.shape,x.shape)
		ax.plot(x[:-1], y)
		ax.axvline(np.mean(all_length), ls = '-', label='{} = {:5.2f}'.format(su,np.mean((all_length))))
	plt.legend(loc=0)
	
	plt.show()
	plt.clf()
	plt.close()


if __name__ =='__main__':
   
    #length_all_subjects()
    length_by_subj()

