import numpy as np
import os
import scipy.io as sio
from sklearn.linear_model import LinearRegression as Linreg
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import Lasso,Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline,Pipeline
from sklearn.model_selection import cross_validate
from sklearn.svm import SVR
from sklearn.model_selection import RandomizedSearchCV, permutation_test_score
import random
random.seed(10)
def GetModel(regname=None,optimisation=False,cv=None):
    if regname.lower()=='ridge':
        model=Ridge()
        if optimisation:
            p_grid = {}
            p_grid['alpha']=np.logspace(-10, 10, 50)
            n_iter_search=50
            Rnd_Srch = RandomizedSearchCV(model,
                                          param_distributions=p_grid,
                                          n_iter=n_iter_search,
                                          cv=cv,iid=True)
            regressor = make_pipeline(StandardScaler(),Rnd_Srch)
        else:
            regressor=Pipeline(steps=[('std',StandardScaler()),('reg',model)])
    if regname.lower()=='lasso':
        model=Lasso()
        if optimisation:
            p_grid = {}
            p_grid['alpha']=np.logspace(-10, 10, 50)
            n_iter_search=50
            Rnd_Srch = RandomizedSearchCV(model,
                                          param_distributions=p_grid,
                                          n_iter=n_iter_search,
                                          cv=cv,iid=True)
            regressor = make_pipeline(StandardScaler(),Rnd_Srch)
        else:
            regressor=Pipeline(steps=[('std',StandardScaler()),('reg',model)])

    if regname.lower()=='linear_model':
        model=Linreg()
        regressor=Pipeline(steps=[('std',StandardScaler()),('reg',model)])


    if regname.lower()=='svr':
        model= SVR(kernel='linear')
        if optimisation:

            p_grid = {}

            p_grid['gamma']= ['scale','auto']
            p_grid['C']= [1, 5,10,15,20,25]
            n_iter_search=10
            Rnd_Srch = RandomizedSearchCV(model, param_distributions=p_grid,
                                           n_iter=n_iter_search, cv=cv,iid=True)
            regressor = make_pipeline(StandardScaler(),Rnd_Srch)
        else:
            regressor=Pipeline(steps=[('std',StandardScaler()),('reg',model)])
    if regname.lower()=='mlp':
        model=MLPRegressor()

        if optimisation:
            p_grid = {}
            p_grid['hidden_layer_sizes']=np.arange(10)+1
            p_grid['activation']=['identity', 'logistic', 'tanh', 'relu']

            n_iter_search=50
            Rnd_Srch = RandomizedSearchCV(model,
                                          param_distributions=p_grid,
                                          n_iter=n_iter_search,
                                          cv=cv,iid=True)
            regressor = make_pipeline(StandardScaler(),Rnd_Srch)
        else:
            regressor=Pipeline(steps=[('std',StandardScaler()),('reg',model)])
    return regressor


def MakeRegression(model=[],X=[],y=[],
                   inner_cv=None,outer_cv=None,
                   feat_selection=False,get_estimator=False,
                   stat=False,nperms=None,
                   njobs=1):
    permutation_scores,pvalue=[],[]
    results=dict()
    if feat_selection==False:
        scores= cross_validate(model, X, y, cv=outer_cv,
                            scoring=('r2', 'neg_mean_squared_error'),
                            return_train_score=False,return_estimator=True)
        results={'r2':scores['test_r2'],'mse':scores['test_neg_mean_squared_error']}
        if get_estimator:
            reslt=dict()
            models=scores['estimator']
            coef=[]

            for model in models:
                if 'randomizedsearchcv' in model.named_steps:
                    coef.append(model.named_steps['randomizedsearchcv'].best_estimator_.coef_)
                else:
                    coef.append(model.named_steps['reg'].coef_)
            results.update({'coeff':coef})
        if stat:
            print('Runing permutations .... ')
            scores,permutation_scores,pvalue=permutation_test_score(model,
                                                                   X,
                                                                   y,
                                                                   scoring='neg_mean_squared_error',
                                                                   cv=outer_cv,
                                                                   n_permutations=nperms,
                                                                   n_jobs=njobs)
            results.update({'scores':scores,'permutation_scores':permutation_scores,'pvalue':pvalue})


    else:
        scores = cross_validate(model, X, y, cv=outer_cv,
                                scoring=('r2', 'neg_mean_squared_error'),
                                return_train_score=False)
        results={'scores':scores}

    return results
