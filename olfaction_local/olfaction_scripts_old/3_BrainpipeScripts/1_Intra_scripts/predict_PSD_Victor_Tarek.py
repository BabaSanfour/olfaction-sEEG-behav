import numpy as np
import scipy.io as sio
import os
import pandas as pd
from regression_utils import MakeRegression,GetModel
from sklearn.model_selection import StratifiedKFold
import random
import json
random.seed(10)
######################## Load PSD and csv files ################################

data_path=os.getcwd()
save_path=os.path.join(data_path,'Results','predict_PSD')
if not os.path.exists(save_path):
    os.makedirs(save_path)
Psd_file=os.path.join(data_path,'PSD_100_ROIs.mat')
PSD_all_bands=sio.loadmat(Psd_file)['PSD']
var_file=os.path.join(data_path,'N_45_full_data_PA_tarek.xls')
df=pd.read_excel(var_file)
X=df[['gender', 'TIme', 'AgeDiag', 'Education', 'Cortico', 'ITMTX', 'RadioDose']]

########################### Test On/off  #######################################

test=True
if test:
    bandes=['Delta']
else:
    bandes=['Delta','Theta','Alpha','Beta','Gamma1','Gamma2','Gamma3','GammaL']

####################### Set regression parameters ##############################

#target_name='Radio'# Set the target name
#target=df[target_name].tolist()
regressor='svr'  # Regressor name
inner_cv=10       # Inner cross validation used to optimise the model's params
outer_cv=10       # Outer cv used to train and test data
optimise=False     # Turn to True if you want to optimise the model
FeatSelect=False  # Set to True if you want to run feature selection
get_estimator=True #set it to True if you want to save all model results
Stat=False        #set stat to true to run permutation tests
n_perms=100    # set the number of permutation to run if stat= True
# [float(i) for i in target]
# print(type(target[0]))
########################### Run Regression #####################################

for bdi,bd in enumerate(bandes):
    print('Runing regression for bande {}'.format(bd))
    y=np.squeeze(PSD_all_bands[:,bdi,:]).T
    model=GetModel(regname=regressor,
                   optimisation=optimise,
                   cv=inner_cv)
    result_all_roi=[]
    for roi in range(y.shape[1]):
        print('Predict ROI nb: {}'.format(str(roi)))
        data=X.values
        target=y[:,roi].astype(float)

        result=MakeRegression(model=model,
                              X=data,
                              y=target,
                              inner_cv=inner_cv,
                              outer_cv=outer_cv,
                              stat=Stat,nperms=n_perms,
                              get_estimator=get_estimator,
                              njobs=-1)
        result_all_roi.append(result)


        print(np.mean(result['r2']))
    save_file=os.path.join(save_path,'Res_100ROI_{b}_{reg}.mat'.format(b=bd,
                                                                reg=regressor))
    print('saving results for band {}'.format(bd))
    sio.savemat(save_file,{'results_ROI':result_all_roi})
