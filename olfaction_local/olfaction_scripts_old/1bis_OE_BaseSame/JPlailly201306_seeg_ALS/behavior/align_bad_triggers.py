# -*- coding: utf-8 -*-

import numpy as np

from itertools import combinations

from scipy import stats

from matplotlib import pyplot


def align_bad_triggers_brut_force(ok, bad, maxiter = 1000):
    
    err = np.inf
    res = None
    for take in combinations(range(bad.size), ok.size):
        #~ print take
        a,b,r,tt,stderr=stats.linregress(ok, bad[list(take)])
        if stderr<err:
            res = list(take)
            err = stderr
            #~ print err
    return np.array(res)
    
def align_bad_triggers_mcmc(ok, bad, maxiter = 1000000, thresh_err = 5e-5):
    mask=  np.ones(bad.size, dtype = bool)
    
    n = ok.size
    err = np.inf
    
    take  = np.arange(bad.size)
    np.random.shuffle(take)
    take = take[:n]
    
    a,b,r,tt,stderr=stats.linregress(ok, bad[take])
    err = stderr
    
    for i in range(maxiter):
        
        r = np.random.randint(n)
        putin = np.arange(bad.size)[~np.in1d(np.arange(bad.size), take)]
        putin = putin[np.random.randint(putin.size)]
        
        putout = take[r]
        take[r] = putin
        bad2 = bad[take]
        bad2.sort()
        a,b,r,tt,stderr=stats.linregress(ok, bad2)
        #~ print stderr
        
        if stderr<err:
            keep = np.random.rand()<.9999
        else:
            keep = np.random.rand()<.0001
        
        if stderr<thresh_err:
            print stderr
            breakpyplot
        
        if  keep:
            err = stderr
            print err
        else:
            take[r] = putout
        
    take.sort()
    return take



def align_bad_triggers_iterative(ok, bad, maxiter = 1000000, thresh_stderr = 1e-3, thresh_slope = 1.e-2):
    
    take = [ 0, 1]
    
    def is_ok():
        if 2<len(take)<=3:
            return stderr<thresh_slope and np.abs(a-1.)<thresh_slope 
        elif len(take)>3:
            #~ return stderr<thresh_slope and np.abs(a-1.)<thresh_slope and r>last_r
            #~ return stderr<thresh_slope and np.abs(a-1.)<thresh_slope and pvalue<last_pvalue
            return stderr<thresh_slope and np.abs(a-1.)<thresh_slope and stderr<last_stderr*2
            #~ return stderr<thresh_slope and np.abs(a-1.)<thresh_slope
        else:
            return np.abs(a-1)<thresh_slope
    
    last_stderr = np.inf
    last_pvalue = 1.
    last_r = 0.
    while True:
        #~ print 
        #~ print take
        #~ print ok
        
        ok2 = ok[:len(take)]
        bad2 = bad[take]
        #~ print np.diff(ok2), np.diff(bad2)
        #~ print np.diff(ok2) - np.diff(bad2)
        
        
        a,b,r,pvalue,stderr = stats.linregress(ok2, bad2)
        
        print 
        print take
        print 'stderr', stderr, last_stderr, stderr<last_stderr
        print 'r', r, last_r, r>last_r
        print 'pvalue', pvalue, last_pvalue, pvalue<last_pvalue
        print is_ok()
        
        #~ err = stderr
        
        if is_ok():
            if len(take)<ok.size:
                if take[-1]<bad.size-1:
                    #~ print 3
                    take.append(take[-1]+1)
                    last_stderr = stderr
                    last_pvalue = pvalue
                    last_r = r
                else:
                    if take[0]+2>=bad.size:
                        #~ print 5
                        return None
                    #~ print 4
                    take = [ take[0]+1, take[0]+2]
                    last_stderr = np.inf
                    last_pvalue = 1.
                    last_r = 0
            else:
                #~ print 6
                break
        else:
            if take[-1]<bad.size-1:
                #~ print 1
                take[-1] += 1
            else:
                #~ print 2
                if take[0]<bad.size-2:
                    take = [ take[0]+1, take[0]+2]
                    last_stderr = np.inf
                    last_pvalue = 1.
                    last_r = 0
                else:
                    #~ print 'pas de solution'
                    return None
        
        #~ print take
        #~ print
    return take



def align_bad_triggers_iterative_v2(ok, bad, maxiter = 1000000, thresh_stderr = 1e-3, thresh_slope = 1.e-2):
    
    take = [ 0, 1]
    #~ last_stderr = np.inf
    #~ last_pvalue = 1.
    #~ last_r = 0.
    
    def is_ok(stderr, a):
        if len(take)==2:
            return np.abs(a-1)<thresh_slope
        else:
            return stderr<thresh_slope and np.abs(a-1.)<thresh_slope
    
    def is_better():
        if len(take)==2:
            return abs(1.-a3)<abs(1.-a2)
        else:
            return pvalue3<pvalue2
            #~ return stderr3<stderr2
            #~ return r3>r2
    
    while True:
        # with last point
        #~ print 
        #~ print take, len(ok), len(bad)
        ok2 = ok[:len(take)]
        bad2 = bad[take]
        a2,b2,r2,pvalue2,stderr2 = stats.linregress(ok2, bad2)
        
        take3  = list(take)
        take3[-1] += 1
        ok3 = ok[:len(take3)]
        bad3 = bad[take3]
        a3,b3,r3,pvalue3,stderr3 = stats.linregress(ok3, bad3)
        
        
        
        #~ print 'a2,b2,r2,pvalue2,stderr2', a2,b2,r2,pvalue2,stderr2
        #~ print take3
        #~ print 'a3,b3,r3,pvalue3,stderr3', a3,b3,r3,pvalue3,stderr3
        #~ print 'is_better()', is_better()
        if is_better():
            take = take3
        else:
            take = take
            #~ print 'is_ok(stderr2, a2)', is_ok(stderr2, a2)
            if is_ok(stderr2, a2):
                if take[-1]<bad.size-1 and len(take)<len(ok):
                    #~ print 'ici add next'
                    take.append(take[-1]+1)
                else:
                    if take[0]+2>=bad.size:
                        #~ print 5
                        return None
                    else:
                        #solution
                        break
            else:
                if take[-1] + 1 == bad.size:
                    return None
                else:
                    take[-1] += 1
        
        if take[-1]+1 == bad.size:
            ok2 = ok[:len(take)]
            bad2 = bad[take]
            a2,b2,r2,pvalue2,stderr2 = stats.linregress(ok2, bad2)            
            if is_ok(stderr2, a2) and len(take) == ok.size:
                break
                #solution
            else:
                #~ print 'lllll', take
                if take[0]<bad.size-3:
                    #new start
                    take = [ take[0]+1, take[0]+2]
                    #~ print 'lllllll', take
                else:
                    #~ print 'pas de solution'
                    return None
    return take




#~ align_bad_triggers = align_bad_triggers_mcmc
#~ align_bad_triggers = align_bad_triggers_iterative
align_bad_triggers = align_bad_triggers_iterative_v2
    
    
    
def test1():
    
    ok = np.array([ 1.,2.,3.,4.])+.3
    bad = np.array([1.,2.,2.2,3.,3.4,4,5.2 ])
    #~ mask = align_bad_triggers_brut_force(ok, bad)
    mask = align_bad_triggers_mcmc(ok, bad)
    
    print mask, ok, bad[mask]
    

def test2():
    
    #~ ok = np.array([ 1.,2.,3.,4.])+.3
    #~ bad = np.array([1.,2.,2.2,3.,3.4,4,5 ])
    ok = np.arange(13)+3
    bad = np.concatenate([np.arange(13), np.arange(5)+.5, np.arange(4)+.4, ])
    bad.sort()
    print bad
    mask = align_bad_triggers_mcmc(ok, bad)
    
    print mask, ok, bad[mask]
    

def test3():
    
    #~ ok = np.array([ 1.,2.,3.,4.])+.3
    #~ bad = np.array([1.,2.,2.2,3.,3.4,4,5 ])
    ok = np.arange(13)+3
    bad = np.concatenate([np.arange(13)+.3, np.arange(5)+.5, np.arange(4)+.4, [-1, 0,.2,.3]])
    bad.sort()
    #~ ok = np.array([   9.67000008,   33.24000168,   51.43999863,   62.88999939,   85.26000214,
                                    #~ 96.01000214,  120.95999908,  143.47999573,  154.71000671,  178.08999634,
                                    #~ 195.44000244, ])
    #~ bad = np.array([  21.65625,   45.25,      63.4375,    74.53125,   74.875,     97.25,     108.,
                                  #~ 132.96875,  155.46875,  165.03125,  165.125,    166.71875,  190.09375,
                                  #~ 207.4375, ])


    #~ ok = np.array([  13.64999962 ,  40.38000107,   70.27999878,   94.84999847,  117.51000214,
                  #~ 144.78999329,  157.78999329,  170.16000366,  182.1000061,   206.27999878,
                  #~ 225.61000061])
    #~ bad = np.array([  24.9375,    39.0625,    39.09375,   39.125,     39.1875,    39.21875,
               #~ 51.65625,   79.65625,   81.5625,   106.15625,  128.8125,   156.09375,
              #~ 163.09375,  163.125,    169.09375,  181.46875,  193.40625,  199.1875,   217.5625,
              #~ 236.90625,  287.25,     325.15625,])

    #~ ok = np.array([  94.30999756,  117.5,         130.88000488,  144.80000305,  167.97000122,
                    #~ 183.80999756,  202.1499939,   219.75999451,  233.02999878,  254.83999634,])
    #~ bad = np.array([  76.96875,  102.6875,   120.4375,   125.90625,  139.28125,  153.1875,   176.375,
                    #~ 192.21875,  210.5625,   228.15625,  241.4375,   252.1875,   263.25,   ])
    
    
    #Date=2013-10-24_Sujet=S02_Exp=E_Run=2_Group=1_Image=M.res
    ok = np.array([6.67999983,  17.92000008,  37.5, 48.66999817,  56.31000137,
            65.34999847,  75.44999695,  87.76999664,  96.26000214, 104.48999786,
            115.51000214, 128.66999817, 135.97999573, 145.3500061,  153.63999939,
            170.75999451, 192.30999756, 212.77999878, 224.47999573, 257.36999512,
            293.97000122, 360.16000366, 435.61999512, 448.39001465, 489.26998901,
            502.10998535, 586.72998047, 594.33001709])
    bad = np.array([3.0625, 26.5625, 37.78125,  57.375,  68.53125,  75.78125,
            75.8125, 76.1875, 85.21875,  95.34375, 107.65625, 116.125, 124.375,
            135.375, 148.5625,  155.84375, 165.21875, 173.53125, 190.625, 201.0625,
            201.09375, 212.1875,  232.65625, 244.34375, 277.25,  313.84375,
            328.90625, 328.96875, 329., 380.03125, 439.40625, 439.4375,  455.5,
            468.28125, 487.65625, 509.15625, 522., 568.15625, 606.625, 613.9375,
            614.21875])
    
    ##Date=2015-10-12_Sujet=S06_Exp=E_Run=1_Group=1_Image=1.res
    #~ ok = np.array([7.53999996,  17.46999931,  39.70000076,  51.83000183,  65.69000244,
            #~ 73.86000061,  82.62999725,  90.98999786, 100.61000061, 108.13999939,
            #~ 117.51000214, 127.40000153, 142.50999451, 149.47000122, 156.05999756,
            #~ 163.3500061,  171.8999939,  179.55999756, 188.13000488, 195.8500061,
            #~ 201.36000061, 209.6000061,  216.88999939, 223.83999634, 232.21000671,
            #~ 246.67999268, 255.46000671, 262.32000732, 268.67999268, 276.17999268,
            #~ 281.80999756, 289.10998535, 297.3500061,  307.70999146, 316.26000977,
            #~ 328.25, 339.44000244, 347.92999268, 355.38000488, 362.63000488,
            #~ 378.04998779, 388.72000122, 396.,403.45001221, 423.57998657,
            #~ 429.79000854, 438.73001099, 445.29998779])
    #~ bad = np.array([15.15625,  117.03125,  117.0625, 155.15625,  257.125, 257.21875,
            #~ 295.5625, 394.15625,  394.25, 519.65625,  519.6875, 543.90625,
            #~ 553.84375,  561.875,  576.0625, 588.21875,  602.0625, 610.25,
            #~ 619., 627.375,  637., 644.53125,  649.65625,  649.75,
            #~ 653.875,  663.78125,  678.875,  685.84375,  691.75, 692.4375,
            #~ 699.71875,  708.28125,  715.9375, 724.5,  732.21875,  737.75,
            #~ 745.96875,  753.28125,  760.21875,  768.59375,  783.0625, 791.84375,
            #~ 798.6875, 805.0625, 812.0625, 812.5625, 818.1875, 825.5,
            #~ 833.71875,  844.09375,  852.65625,  864.625,  875.8125, 884.3125,
            #~ 891.75, 897.75, 897.78125, 897.84375,  899., 914.4375,
            #~ 925.09375,  932.375,  939.84375,  946., 959.96875,  966.1875,
            #~ 975.125,  981.6875,  1016.4375, ])
    
    #~ ##Date=2015-10-13_Sujet=S06_Exp=E_Run=2_Group=1_Image=1.res
    #~ ok = np.array([5.53999996,  16.97999954,  30.43000031,  44.04999924,  57.84000015,
        #~ 67.26000214,  78.69000244,  87.88999939,  98.55000305, 107.01000214,
        #~ 115.23999786, 122.80000305, 128.8999939,  135.99000549, 142.08000183,
        #~ 150.86000061, 158.16999817, 167.02000427, 173.77000427, 180.82000732,
        #~ 189.38999939, 196.11999512, 201.21000671, 207.8999939,  214.8500061,
        #~ 221.41999817, 233.25, 238.00999451, 244.83000183, 250.13999939,
        #~ 257.07000732, 262.89001465, 269.16000366, 276.3999939,  283.14001465,
        #~ 289.79998779, 296.79998779, 303.20999146, 310.05999756, 317.95999146,
        #~ 326.69000244])
    #~ bad = np.array([91.8125, 91.9375,  135.1875,  135.21875, 214.09375, 214.125, 214.21875,
        #~ 257.5625,  336.0625,  336.09375, 373.5625,  379.125, 385.,  398.4375,
        #~ 412.0625,  425.84375, 435.28125, 446.6875,  455.90625, 466.5625,
        #~ 473.46875, 473.59375, 475.03125, 483.25,  490.8125,  496.90625, 504.,
        #~ 510.09375, 518.875, 526.1875,  535.03125, 541.78125, 548.84375,
        #~ 557.40625, 564.125, 569.21875, 575.90625, 582.875, 589.4375,  601.25,
        #~ 606.03125, 612.84375, 618.15625, 625.09375, 630.90625, 637.1875,
        #~ 644.40625, 651.15625, 655.90625, 657.8125,  664.8125,  671.21875,
        #~ 678.0625,  685.96875, 694.71875, 742.6875,  786.25])
        
    #~ Date=2015-10-14_Sujet=S06_Exp=R_Run=1_Group=1_Rand=1_Ordre=.res
    #pas de solution
    #~ ok = np.array([14.19999981,  30.18000031,  58.38999939,  88.15000153, 113.19000244,
            #~ 141.41000366, 180.80999756, 201.02999878, 220.00999451, 235.82000732,
            #~ 265.66000366, 297.86999512, 321.45001221, 338.04998779, 357.80999756])
    #~ bad = np.array([39.46875, 125.1875,  125.21875, 166.,  240.75,  240.8125,
            #~ 281.46875, 369.40625, 369.4375,  369.5625,  494.8125,  494.90625,
            #~ 534.1875,  594.40625, 594.4375,  635.46875, 648.75])
    
    #~ Date=2015-10-14_Sujet=S06_Exp=R_Run=3_Group=1_Rand=1_Ordre=.res
    #~ ok = np.array([15.94999981,  30.07999992,  51.15000153,  77.59999847,  89.15000153,
            #~ 105.98000336, 124.12000275, 140.63999939, 151.44000244, 174.61000061,
            #~ 197.33999634, 210.97000122, 236.63999939, 250.02000427,])
    #~ bad = np.array([46.875,  46.9375, 90.4375,  102.875, 117.,  138.0625,
        #~ 164.53125, 176.03125, 176.0625,  192.90625, 211.03125, 221.03125,
        #~ 227.5625,  238.34375, 261.53125, 284.25,  297.875, 312.34375, 312.375,
        #~ 312.4375,  323.5625,  336.9375,  356.875, ])
    
    
    #~ print bad
    #~ mask = align_bad_triggers_mcmc(ok, bad)
    #~ mask = align_bad_triggers_iterative(ok, bad)
    mask = align_bad_triggers_iterative_v2(ok, bad)
    print mask


    fig, ax = pyplot.subplots()
    for i in ok:
        ax.axvline(i, color = 'b')
    for i in bad:
        #~ if i<bad[mask][0] or i>bad[mask][-1]: continue
        ax.axhline(i, color = 'c')

    
    if mask is None:
        print 'pas de solution'
    else:
        #~ print ok - bad[mask]
        for i in bad[mask]:
            ax.axhline(i, color = 'b')
        
        a,b,r,tt,stderr=stats.linregress(ok, bad[mask])
        
        ax.plot(ok, ok*a+b, color = 'r')
    
        
    #~ print np.diff(ok)
    #~ print np.diff(bad[mask])
        #~ print np.diff(ok) - np.diff(bad[mask])
    ax.set_aspect('equal')
    
    pyplot.show()
    
    

if __name__ == '__main__':
    #~ test1()
    #~ test2()
    test3()

    
