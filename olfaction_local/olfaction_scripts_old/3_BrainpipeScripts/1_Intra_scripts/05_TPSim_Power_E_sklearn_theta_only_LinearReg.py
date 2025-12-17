from os.path import join
from itertools import product
import numpy as np
import scipy.io as sio

from brainpipe.system import study
from utils import rename_elecs
from regression_utils import MakeRegression, GetModel
from sklearn.model_selection import StratifiedKFold
import random
random.seed(10)

#################### Set Parameters for the data #################################"
exp, freq = 'Enc', 'theta'
meth, conds = 'btw', ['low','mid','high']
subjects = ['CHAF','VACJ','PIRJ','SEMC','FERJ','LEFC']
rois_to_keep = ['aHC','HC','MFG','ACC','IFG','Amg','pPirT','PHG','Ins_olf','OFC_olf','SFG']

#################### Data and save PATHS #################################"
st = study('Olfacto')
path_tps = join(st.path, 'feature/TPSim_3groups_'+exp+'/')
tps_form = join(path_tps, 'TPS_'+meth+'/TPS_pears_{}_{}_'+meth+'_{}.npz')
path_save = join(st.path, 'classified/TPSim_LinReg_{}_3groups_btw_theta/')
name_classif = join(path_save, '{}_{}_'+meth+'_3gr_{}.npz') #su, model, freq

#################### ML parameters #################################"
regressor = 'svr' #svr, lasso, ridge
inner_cv=10       # Inner cross validation used to optimise the model's params
outer_cv=10       # Outer cv used to train and test data
optimise=False     # Turn to True if you want to optimise the model
FeatSelect=False  # Set to True if you want to run feature selection
get_estimator=True #set it to True if you want to save all model results
Stat = True #set stat to true to run permutation tests
n_perms = 1000

#################### Run Regression #################################"

for su in subjects:
    mat0 = np.load(tps_form.format(su,conds[0],freq),allow_pickle=True)
    labels = mat0['label']
    x,y,z = mat0['xyz'][:,0], mat0['xyz'][:,1], mat0['xyz'][:,2]
    new_labels = rename_elecs(labels,x,y,z)
    idx_rois = np.isin(new_labels, rois_to_keep)
    nelecs = mat0['tps'][idx_rois].shape[0]
    
    tps_list, perf_code = [], []
    for i,cond in enumerate(conds):
        data = np.load(tps_form.format(su,cond,freq))['tps'][idx_rois]
        tps_list.append(data)
        perf_code.extend([i+1]*data.shape[-1])
    print (su,mat0.files, 'TPS shape: ', [tps.shape for tps in tps_list],len(perf_code))
    
    #=========================== Create dict for all results =================================    
    kwargs = {}
    kwargs['names'], kwargs['channels'] = new_labels[idx_rois], mat0['channel'][idx_rois]
    kwargs['xyz'] = mat0['xyz'][idx_rois]
    kwargs['model'] = regressor

    # =========================== Select Power for 1 elec 1 freq =================================                 
    model = GetModel(regname=regressor, optimisation=optimise,cv=inner_cv)
    regress_results = []
    for elec_num in range(nelecs):
        print('--» processing',su, 'elec', elec_num,'/',nelecs)
        tps_data_elec = [tps[elec_num][:,np.newaxis] for tps in tps_list]

        # create a data matrix, concatenate along the trial dimension
        x = np.concatenate(tps_data_elec, axis=0)
        y = np.array(perf_code)
        print ('Size of the concatenated data: ', x.shape)
        results = MakeRegression(model=model,X=x,y=y,inner_cv=inner_cv,
                                outer_cv=outer_cv, stat=Stat, nperms= n_perms,
                                get_estimator=get_estimator,njobs=-1)
        regress_results.append(results)
        print('R2 mean', np.mean(results['r2']))
        print('score mean', np.mean(results['scores']))
        print('permut mean', np.mean(results['permutation_scores']))
    
    kwargs['results'] = regress_results
    #Save plots
    np.savez(name_classif.format(exp,su,regressor,freq), **kwargs)
    del x, y, tps_data_elec, results, regress_results
    del tps_list