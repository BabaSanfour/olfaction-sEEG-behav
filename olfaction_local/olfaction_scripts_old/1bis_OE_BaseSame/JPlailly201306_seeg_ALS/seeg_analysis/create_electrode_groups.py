# -*- coding: utf-8 -*-


import sys, os
sys.path.append('../behavior') 

from connection import *

from collections import OrderedDict

from PyQt4 import QtCore, QtGui
import pyqtgraph as pg

import pyplotbrain as ppb


def create_electrode_groups(force = True):
    for subject in session.query(Subject):
        bl = subject.blocks[0]
        rcgs = bl.recordingchannelgroups
        rcs = rcgs[0].recordingchannels.order_by('index').all()
        
        # if alreay done
        if rcgs.count()>1:
            if not force:
                continue
            #delete
            for rcg in rcgs:
                if rcg.name == 'all channels': continue
                session.delete(rcg)
            session.commit()
        
        letters, nums = [ ], [ ]
        for rc in rcs:
            if ' ' in rc.name:
                letter, num = rc.name.split(' ')
            else:
                letter, num = None, None
            letters.append(letter)
            nums.append(num)
        letters = np.array(letters, dtype = str)
        
        for letter in np.unique(letters):
            if letter == 'None': continue
            rcg = RecordingChannelGroup(name = letter)
            rcgs.append(rcg)
            
            for i in np.where(letters == letter)[0]:
                rcg.recordingchannels.append(rcs[i])
        session.commit()




def test_plot_bipolar():
    subject = session.query(Subject)[0]
    bl = subject.blocks[0]
    rcgs = bl.recordingchannelgroups
    
    rcg = rcgs.filter_by(name = 'B')[0]
    rcs = rcg.recordingchannels.order_by('`index`').all()
    
    sigs = np.vstack([rc.analogsignals[0].signal.magnitude for rc in rcs])
    
    from matplotlib import pyplot
    from OpenElectrophy.timefrequency import TimeFreq
    n = sigs.shape[0]
    
    
    
    pts = 50000
    sr = rc.analogsignals[0].sampling_rate.rescale('Hz').magnitude
    times = np.arange(pts)/sr
    
    
    f_start = 1.
    f_stop = 120.
    clim = [0, 5]
    
    fig1, axs = pyplot.subplots(nrows = n,ncols = 2, sharex = True, sharey = True)
    fig2, ax2s = pyplot.subplots(nrows = n,ncols = 2, sharex = True, sharey = True)
    for i in range(n):
        axs[i, 0].plot(times, sigs[i, :pts])
        
        neo_sig = rcs[i].analogsignals[0].to_neo()
        neo_sig = neo_sig[:pts]
        #~ print neo_sig.sampling_rate
        tfr = TimeFreq(neo_sig, f_start = f_start, f_stop = f_stop, sampling_rate = sr*pq.Hz)
        tfr.plot(ax = ax2s[i, 0], clim =clim, colorbar = False)
        
        if i<sigs.shape[0]-1:
            axs[i, 1].plot(times, sigs[i, :pts] - sigs[i-1, :pts])
            neo_sig2 = rcs[i+1].analogsignals[0].to_neo()
            neo_sig2 = neo_sig2[:pts]
            tfr = TimeFreq(neo_sig-neo_sig2, f_start = f_start, f_stop = f_stop, sampling_rate = sr*pq.Hz)
            tfr.plot(ax = ax2s[i, 1], clim = clim, colorbar = False)
            
            
            
    axs[n-1, 1].plot(times, np.zeros(pts))
    
    pyplot.show()
    

    

#~ def create_bipolar(force = True):
    #~ for subject in session.query(Subject):
        #~ bl = subject.blocks[0]
        #~ rcgs = bl.recordingchannelgroups
        
        #~ # if alreay done
        #~ if rcgs.filter(RecordingChannelGroup.name.like('Bipolar')).count():
            #~ # TODO
            #~ if not force : continue
            #~ #delete
            #~ print rcgs.filter_by(RecordingChannelGroup.name.like('Bipolar')).all()
        
        #~ for rcg in rcgs:
            #~ if rcg.name == 'all channels': continue


    
if __name__ =='__main__':
    create_electrode_groups()

    
    #test_plot_bipolar()
    
