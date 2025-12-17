from mne import *
import numpy as np
from os import path
from brainpipe.system import study
from brainpipe.feature import power, amplitude, sigfilt


st = study('Olfacto')
pathdata = path.join (st.path, 'database/TS_E_all_cond_by_block_trigs_filter1_500art/')
save2path = path.join(st.path, 'features/TS_E_all_cond_by_block_trigs_filter1_500art/')

    
# Files to load
subj = ['CHAF', 'VACJ', 'SEMC', 'FERJ']
trigger =  ['odor', 'rest',]
fnames = ['delta', 'theta', 'alpha', 'beta', 'gamma30-60', 'gamma60-120']
freq = {	
			'delta' : [2, 4], 
			'theta' : [4, 8], 
			'alpha' : [8, 13], 
			'beta' : [13, 30], 
			'gamma30-60' : [30, 60], 
			'gamma60-120' : [60, 120],
		}

for su in subj:
    for trigg in trigger:
    	file = su+'_'+trigg+'_E1E2_concat_allfilter1_bipo.npz'
    	for f in fnames:
    		# Define power settings :
	        kwargs = {} # Define an empty dictionnary to save all power parameters
	        print (f)
	        kwargs['f'] = freq[f] # Frequency vector
	        fname = f
	        #kwargs['baseline'] = (10, 1536) # Where your baseline (start, end) (IN SAMPLE) rest period ~ 3s
	        #kwargs['norm'] = 0 # Type of normalisation (see help on power for more details)
	        kwargs['width'], kwargs['step'] = 100, 50 # take power in 100 samples windows width every 50 samples
	        # Load file :
	        loadname = path.join(st.path, 'database/TS_E_all_cond_by_block_trigs_filter1_500art/', file)
	        print('-> Compute power on: '+loadname)
	        mat = np.load(loadname)
	        x, sf = mat['x'], int(mat['sf'])
	        n_elec, n_pts, n_trials = x.shape
	        # Define a power object :
	        powObj = power(sf, n_pts, **kwargs)
	        kwargs['xpow'],  kwargs['xpow_pval']= powObj.get(x, n_jobs=-1)
	        # Finally save it :
	        kwargs['fname'] = fname
	        savename = loadname.replace('.npz', '_power_'+fname+'.npz').replace('database', 'feature')
	        np.savez(savename, **kwargs)
	        del kwargs, powObj, mat, x, sf, n_elec, n_trials, n_pts
