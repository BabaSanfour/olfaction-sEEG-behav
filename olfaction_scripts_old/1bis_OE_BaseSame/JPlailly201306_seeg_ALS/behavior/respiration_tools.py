# -*- coding: utf-8 -*-

from connection import *

from collections import OrderedDict

# TODO change when OE3 is out
from simplerespiration import SimpleRespirationCycleExtraction




def compute_cycle_for_run(resp):
        detector = SimpleRespirationCycleExtraction()
        
        def t(self,):
            return self.t_start + np.arange(self.signal.shape[0])/self.sampling_rate
        resp2 = type('OldResp', (object,), {'t':t })()
        resp2.signal = resp.signal.magnitude
        resp2.sampling_rate = resp.sampling_rate.magnitude
        resp2.t_start = resp.t_start.magnitude
        
        
        #~ manual_baseline = -1.5
        manual_baseline = 0.
        
        cycle_times = detector.compute(resp2,
        
                                    # preprocessing
                                    inspiration_sign = '-',
                                    high_pass_filter = 12.,
                                    constrain_frequency = 1.,
                                    median_windows_filter = 0.07,
                                    
                                    # baseline
                                    baseline_with_average = True,
                                    #~ manual_baseline =manual_baseline,
                                    
                                    # clean
                                    eliminate_time_shortest_ratio = 5,
                                    eliminate_amplitude_shortest_ratio = 5,
                                    eliminate_mode = 'OR', # 'AND'*
                                    )
        #~ resp.cycle_times = cycle_times*pq.s
        #~ resp.update(session = session)
        return cycle_times*pq.s


def get_time_first_inspi_or_trig(trial, resp, max_percent = .6):
    if resp.cycle_times is None:
        print 'compute resp'
        resp.cycle_times = compute_cycle_for_run(resp)
        resp.save(session)
    
    cycles = resp.cycle_times.magnitude
    
    insp = cycles[:-1,0]
    insp_1 = cycles[1:,0]
    expi = cycles[:-1,1]
    
    t = trial.triggered_odor_time
    #between expi and inspi
    ind, = where((t>expi) & (t<=insp_1))
    if ind.size == 1:
        #~ print 'avant inspi'
        return insp_1[ind[0]]
    
    ind, = where((t>insp) & (t<=expi))
    #~ print t
    #~ print ind
    if ind.size == 1:
        percent = (t-insp[ind[0]])/(expi[ind[0]]-insp[ind[0]])
        
        #~ print 'apres inspi', percent, 
        if percent <max_percent:
            #~ print 'percent OK'
            return t
        else:
            if ind[0]+1<insp.size:
                #~ print 'percent false next one'
                return insp[ind[0]+1]
            else:
                print 'Bug pas de cycle pour trial', trial.id, trial.index, trial.run.subject.name, trial.run.index
                return None

def get_index_cycle(t, resp, max_percent=0.6):

    cycles = resp.cycle_times.magnitude
    insp = cycles[:-1,0]
    insp_1 = cycles[1:,0]
    expi = cycles[:-1,1]
    
    ind, = where((t>expi) & (t<=insp_1))
    if ind.size == 1:
        #between expi and inspi
        ind_cycle = ind[0]+1
    else:
        #between  inspi and expi
        ind, = where((t>insp) & (t<=expi))
        if ind.size == 1:
            percent = (t-insp[ind[0]])/(expi[ind[0]]-insp[ind[0]])
            if percent <max_percent:
                ind_cycle = ind[0]
            else:
                if ind[0]+1<insp.size:
                    ind_cycle = ind[0]+1
                else:
                    ind_cycle = None
        else:
            ind_cycle = None
    
    if ind_cycle == insp.size:
        ind_cycle = None
    
    return ind_cycle
    
    

if __name__ =='__main__':
    #~ test1()
    pass
    