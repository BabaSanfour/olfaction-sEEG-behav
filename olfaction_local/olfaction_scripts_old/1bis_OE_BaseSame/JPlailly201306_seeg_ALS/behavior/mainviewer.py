# -*- coding: utf-8 -*-
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg

from connection import *

from OpenElectrophy.gui.viewers.multiviewer import MultiViewer,SubViewer, SignalViewer,  EpochViewer, EventViewer, TimeFreqViewer, EventList
from OpenElectrophy.gui.viewers.tools import ViewerBase, find_best_start_stop
from OpenElectrophy.gui.guiutil import SpinAndSliderWidget

"""
TODO:
  * color
  * filters

"""


class MainViewer(MultiViewer) :
    def __init__(self, run, parent = None,):
        MultiViewer.__init__(self, parent)
        #~ self.setWindowTitle(u'{}'.format(run_id))
        
        #read data
        self.run = run
        seg = self.neo_seg = run.segments[0].to_neo(cascade = True)
        
        
        self.xsize_changer = SpinAndSliderWidget()
        self.timeseeker.toolbar.addWidget(QtGui.QLabel('xsize:'))
        self.timeseeker.toolbar.addWidget(self.xsize_changer)
        self.xsize_changer.sigChanged.connect(self.change_xsize)
        
        #Signals
        self.add_analogsignals( analogsignals = seg.analogsignals, name = 'Full band signal' )

        #Artefact
        
        artefacts = self.run.segments[0].epocharrays.filter(EpochArray.name.like('Artifact %')).all()
        neo_artefact = [ ea.to_neo(cascade = False) for ea in artefacts]
        self.add_epochs( epocharrays = neo_artefact, name = 'Artifact')
        
        dock1 = self.subviewers[0].dock
        dock2 = self.subviewers[1].dock
        self.tabifyDockWidget(dock1, dock2)
        dock2.setVisible(False)
        
        #resp
        resp = self.run.respirationsignals[0]
        neo_resp = neo.AnalogSignal(resp.signal, t_start = resp.t_start, sampling_rate = resp.sampling_rate, name = 'Resp', color = 'w')
        self.add_analogsignals( analogsignals = [neo_resp], name = 'Resp' )
        epoch_inspi = neo.EpochArray(times = resp.cycle_times[:-1, 0], durations = resp.cycle_times[:-1, 1]-resp.cycle_times[:-1, 0], color = 'g', name = 'inspi')
        epoch_expi = neo.EpochArray(times = resp.cycle_times[:-1, 1], durations = resp.cycle_times[1:, 0]-resp.cycle_times[:-1, 1], color = 'r', name = 'expi')
        self.add_epochs( epocharrays = [epoch_expi, epoch_inspi, ], name = 'cycles')
        
        #TFR
        self.add_timefreqs( analogsignals = seg.analogsignals, name = 'tfr' )
        
        #epoch in DB
        epochs_without_artefact = self.run.segments[0].epocharrays.filter(~EpochArray.name.like('Artifact %')).all()
        epochs_without_artefact = [ ea.to_neo(cascade = False) for ea in epochs_without_artefact]
        epocharrays_behavior = epochs_without_artefact
        self.add_epochs( epocharrays = epocharrays_behavior, name = 'behavior1')
        
        
        # Trials events
        trials = self.run.trials.order_by(asc(text('`Trial`.`index`'))).all()
        labels = [ u'{} odor{} {} {}'.format(t.index, t.odor_num, t.score_recognition, t.score_episodic_strict) for t in trials]
        
        trial_time = neo.EventArray(times = [t.time for t in trials]*pq.s, name = 'trial_time', labels =  labels, color = 'w')
        if run.exp == 'E':
            wanted_odor_time = neo.EventArray(times = []*pq.s, name = 'wanted_odor_time', labels =  [], color = 'c')
        else:
            wanted_odor_time = neo.EventArray(times = [t.wanted_odor_time for t in trials]*pq.s, name = 'wanted_odor_time', labels =  labels, color = 'c')
        triggered_odor_time = neo.EventArray(times = [t.triggered_odor_time for t in trials]*pq.s, name = 'triggered_odor_time', labels =  labels, color = 'm')
        first_inspiration_time = neo.EventArray(times = [t.first_inspiration_time for t in trials]*pq.s, name = 'first_inspiration_time', labels =  labels, color = 'g')
        
        events = [trial_time,wanted_odor_time, triggered_odor_time, first_inspiration_time]
        
        
        events_db = self.run.segments[0].eventarrays.all()
        events_db = [ ev.to_neo(cascade = False) for ev in events_db]
        
        self.add_events( eventarrays = events+events_db, name = 'behavior2')
        # 
        
        
        self.add_eventlist( eventarrays = events + epocharrays_behavior+events_db)
        
        
        
        self.timeseeker.set_start_stop(*find_best_start_stop(segment =seg))
        self.timeseeker.seek(0.)
        
        self.change_xsize(xsize = 20., force = True)

    def change_xsize(self, xsize = None, force = False):
        if xsize is None:
            xsize = self.xsize_changer.value()
        else:
            self.xsize_changer.setValue(xsize)
        for subviewer in self.subviewers:
            if subviewer.dock.isVisible() or force:
                subviewer.viewer.xsize = xsize

    




if __name__ == '__main__':
    #~ run = session.query(Run).get(368)
    #~ run = session.query(Run).get(369)
    run = session.query(Run).get(428)
    
    app = pg.mkQApp()
    w = MainViewer(run = run)
    w.show()
    app.exec_()
