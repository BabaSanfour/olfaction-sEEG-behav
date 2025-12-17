# -*- coding: utf-8 -*-

from connection import *
from collections import OrderedDict
from pandas import DataFrame, Series, concat, ExcelWriter,merge, MultiIndex

import numpy as np

savepath = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/respi/cycles_by_cond/'


for exp in ['E']:
    for subject in session.query(Subject).order_by('Subject.`group`'): #filter_by(name='LEFC'):
        for i,run in enumerate(subject.runs.filter_by(exp=exp).order_by('index')): #filter_by(index=2)
            idx_click, idx_odor0, idx_odor1 = [], [], []
            print subject.name, exp, run.index
            for trial in run.trials:
                t0 = trial.time
                t1 = trial.triggered_odor_time
                t2 = trial.triggered_odor_time + 7.
                resp = trial.run.respirationsignals[0]
                    
                all_cycles = resp.cycle_times.magnitude
                n_cycles = all_cycles.shape[0]-1
                insp = all_cycles[:-1,0]
                expi = all_cycles[:-1,1]
                
                click_cycles, = np.where((insp>t0-1)&(insp<t1))
                odor_cycles, = np.where((insp>t1)&(insp<t2)) #cycles while odor is sent
                if len(odor_cycles) > 0 and len(click_cycles) > 0:
                    odor0, odor1 = odor_cycles[0], odor_cycles[1] if len(odor_cycles)>1 else None
                    idx_odor0.append(odor0)
                    idx_odor1.append(odor1) if odor1 is not None else idx_odor1
                    idx_click.extend(list(click_cycles)) if odor0 not in click_cycles else idx_click
                if len(odor_cycles) > 0 and len(click_cycles) < 0:
                    odor0, odor1 = odor_cycles[0], odor_cycles[1] if len(odor_cycles)>1 else None
                    idx_odor0.append(odor0)
                    idx_odor1.append(odor1) if odor1 is not None else idx_odor1
                if len(odor_cycles) < 0 and len(click_cycles) > 0:
                    idx_click.extend(list(click_cycles))
            #remove duplicates between odor and click cycles (keep in odor)
            print len(idx_click),len(idx_odor0),len(idx_odor1)
            idx_click = [idx for idx in idx_click if idx not in idx_odor1]
            print len(idx_click),len(idx_odor0),len(idx_odor1)

            #create boolean vectors and save them
            bool_click, bool_od0 = np.ones(n_cycles)*0,np.ones(n_cycles)*0
            bool_od1, bool_no_od = np.ones(n_cycles)*0,np.ones(n_cycles)
            bool_od_all = np.ones(n_cycles)*0
            print len(bool_click), len(bool_no_od)
            bool_click[idx_click]=1
            bool_od0[idx_odor0]=1
            bool_od1[idx_odor1]=1
            bool_od_all[idx_odor0+idx_odor1]=1
            bool_no_od[idx_click+idx_odor1+idx_odor0]=0 
            #print sum(bool_click),sum(bool_no_od),sum(bool_od0)+sum(bool_od1)
            #print sum([bool_click+bool_no_od+bool_od0+bool_od1]), n_cycles
            print 'no odors', sum(bool_no_od)
            print 'all odors', sum(bool_od_all)
            
            np.save(savepath+subject.name+'_'+exp+str(run.index)+'_click.npy',bool_click)
            np.save(savepath+subject.name+'_'+exp+str(run.index)+'_no_odor.npy',bool_no_od)
            np.save(savepath+subject.name+'_'+exp+str(run.index)+'_odor0.npy',bool_od0)
            np.save(savepath+subject.name+'_'+exp+str(run.index)+'_odor1.npy',bool_od1)
            np.save(savepath+subject.name+'_'+exp+str(run.index)+'_odorall.npy',bool_od_all)
