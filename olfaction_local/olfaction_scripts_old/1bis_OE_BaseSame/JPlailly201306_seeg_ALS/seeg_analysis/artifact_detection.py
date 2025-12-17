# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior')
from connection import *

import quantities as pq
import numpy as np
import scipy.signal
from OpenElectrophy.timefrequency import TimeFreq

from insertdata import signalToTrig


from matplotlib import pyplot
import time
import neo

import mycolors

from sqlalchemy import inspect



def detect_artifact_on_one_sig(neosig,f1=60, f2=250.,  method = 'timefreq', thresh_factor = 40):
    sr =neosig.sampling_rate.rescale('Hz').magnitude # sampling rate original
    sig = neosig.magnitude
    times = neosig.times.magnitude

    if method == 'timefreq':
        # method1 : tfr
        tfr = TimeFreq(neosig,
                                        t_start = neosig.t_start,
                                        t_stop =neosig.t_stop,
                                        f_start = f1 * pq.Hz,
                                        f_stop = f2* pq.Hz,
                                        deltafreq = (f2-f1)/20.* pq.Hz,
                                        f0 = 2.5,
                                        sampling_rate =  sr * pq.Hz,
                                        #optimize_fft = True,
                                        )
        ampl = np.abs(tfr.map).mean(axis=1)
        med = np.median(ampl)
        mad = np.median(np.abs(ampl-med))
        thresh = med + thresh_factor*mad
        #~ print 'thresh', thresh

        ind1,  ind2 = signalToTrig(ampl, thresh, front = True, clean = True)
        artefact_times = tfr.times.magnitude[ind1]
        artefact_durations = tfr.times.magnitude[ind2] - tfr.times.magnitude[ind1]


    elif method == 'hilbert':
        #~ n_init = int(ana.size*sr/ana.sampling_rate) # normal number of points for signal, after subsampling
        #~ n=2**np.ceil(np.log(n_init)/np.log(2))

        # method 2 : envelope sur filtre passe haut
        n_decimate = int(np.floor(sr/f2/2.))
        print 'n_decimate', n_decimate
        sig2 = scipy.signal.decimate(sig, n_decimate)
        sr2 = sr/n_decimate # sub sample rate
        print 'sr2', sr2
        times2 = np.arange(sig2.size)/sr2

        Wn = [f1/(sr2/2.), f2/(sr2/2.) ]
        b, a = scipy.signal.iirfilter(N=3, Wn=Wn, btype = 'bandpass', analog = 0, ftype = 'butter', output = 'ba')
        sigf = scipy.signal.filtfilt(b, a, sig2)
        #~ print 'yep'
        n_init = sig2.size
        n=int(2**np.ceil(np.log(n_init)/np.log(2)))
        #~ print n_init, n
        ampl2 = np.abs(scipy.signal.hilbert(sigf, N = n))# estimation with hilbert
        #~ printRecordingChannel ampl2.size
        ampl2 = ampl2[:n_init]
        #~ print 'yop'
        med2 = np.median(ampl2)
        mad2 = np.median(np.abs(ampl2-med2))
        thresh2 = med2 + thresh_factor*mad2

        ind1,  ind2 = signalToTrig(ampl2, thresh2, front = True, clean = True)
        artefact_times = times2[ind1]
        artefact_durations = times2[ind2] - times2[ind1]

    return artefact_times, artefact_durations



def test_detect_artifact_on_one_sig():
    #~ run = session.query(Run)[0]
    run = session.query(Run).get(206)
    #~ anasig = run.segments[0].analogsignals[0]
    #~ anasig = run.segments[0].analogsignals.filter_by(channel_index = 0)[0]
    anasig = run.segments[0].analogsignals.filter_by(name = 'B 1')[0]
    print anasig
    neosig = anasig.to_neo()
    print neosig.sampling_rate

    t1 = time.time()
    artefact_times1, artefact_durations1 = detect_artifact_on_one_sig(neosig, method = 'timefreq', f1=30, f2=250.,  thresh_factor = 40.)
    t2 = time.time()
    print 'timefreq', t2-t1, artefact_times1.size

    t1 = time.time()
    artefact_times2, artefact_durations2 = detect_artifact_on_one_sig(neosig, method = 'hilbert', f1=30, f2=250.,  thresh_factor = 40.)
    t2 = time.time()
    print 'hilbert', t2-t1, artefact_times2.size


    fig, axs = pyplot.subplots(nrows = 2, sharex = True)
    axs[1].plot(neosig.times.magnitude, neosig.magnitude, color = 'b')

    for t, d in zip(artefact_times1, artefact_durations1):
        axs[1].axvspan(t, t+d, color = 'm', alpha = .2)

    for t, d in zip(artefact_times2, artefact_durations2):
        axs[1].axvspan(t, t+d, color = 'g', alpha = .2)


    ##DEBUG hilbert run 173
    #~ f1=60
    #~ f2=250.
    #~ thresh_factor = 50.
    #~ sr =neosig.sampling_rate.rescale('Hz').magnitude # sampling rate original
    #~ sig = neosig.magnitude
    #~ times = neosig.times.magnitude
    #~ n_decimate = int(np.floor(sr/f2/2.))
    #~ sig2 = scipy.signal.decimate(sig, n_decimate)
    #~ sr2 = sr/n_decimate # sub sample rate
    #~ times2 = np.arange(sig2.size)/sr2

    #~ Wn = [f1/(sr2/2.), f2/(sr2/2.) ]
    #~ b, a = scipy.signal.iirfilter(N=3, Wn=Wn, btype = 'bandpass', analog = 0, ftype = 'butter', output = 'ba')
    #~ sigf = scipy.signal.filtfilt(b, a, sig2)
    #~ n_init = sig2.size
    #~ n=2**np.ceil(np.log(n_init)/np.log(2))
    #~ ampl2 = np.abs(scipy.signal.hilbert(sigf, N = n))# estimation with hilbert
    #~ ampl2 = ampl2[:n_init]
    #~ med2 = np.median(ampl2)
    #~ mad2 = np.median(np.abs(ampl2-med2))
    #~ thresh2 = med2 + thresh_factor*mad2

    #~ axs[0].plot(times2, ampl2)
    #~ axs[0].axhline(thresh2)


    pyplot.show()



def compute_and_store_artifact(anasig, dbinfo = None):
    name =  'Artifact {} {}'.format(anasig.name, anasig.id)
    session = inspect(anasig).session

    if dbinfo is not None:
        EpochArray_ = dbinfo.classes_by_name['EpochArray']
        ep = session.query(EpochArray_).filter_by(name = name).first()

    ep = session.query(EpochArray).filter_by(name = name).first()
    if ep is not None:
        session.delete(ep)
        session.commit()


    neosig = anasig.to_neo()
    artefact_times, artefact_durations = detect_artifact_on_one_sig(neosig,f1=30, f2=250.,  method = 'hilbert', thresh_factor = 20.)
    ep = EpochArray(name = name, times = artefact_times*pq.s, durations = artefact_durations* pq.s)
    anasig.segment.epocharrays.append(ep)
    session.commit()

    return ep

def test_compute_and_store_artifact():
    run = session.query(Run)[0]
    print run.subject

    for anasig in  run.segments[0].analogsignals:
        print anasig.name, anasig.signal.size
        compute_and_store_artifact(anasig)


def is_artifact_in_window(anasig, t1,t2, dbinfo = None):
    name =  'Artifact {} {}'.format(anasig.name, anasig.id)
    session = inspect(anasig).session
    ep = session.query(EpochArray).filter_by(name = name).first()
    if ep is None:
        print 'Compute artifact for', anasig.name, anasig.id,
        ct1 = time.time()
        ep = compute_and_store_artifact(anasig, dbinfo = dbinfo)
        ct2 = time.time()
        print '   during', ct2-ct1, 's.'
    #~ print ep
    starts, stops = ep.times.magnitude, ep.times.magnitude+ep.durations.magnitude
    #inside
    ok1  = np.any(((starts<=t1)&(stops>t1)) | ((starts<=t2)&(stops>t2)))
    # or contain
    ok2 = np.any((starts>=t1 ) & (stops<t2 ))

    return ok1 or ok2

def return_artifact_timings(anasig, dbinfo = None):

    name =  'Artifact {} {}'.format(anasig.name, anasig.id)
    session = inspect(anasig).session
    ep = session.query(EpochArray).filter_by(name = name).first()
    if ep is None:
        print 'Compute artifact for', anasig.name, anasig.id,
        ct1 = time.time()
        ep = compute_and_store_artifact(anasig, dbinfo= dbinfo)
        ct2 = time.time()
        print '   during', ct2-ct1, 's.'
    #~ print ep
    starts, stops = ep.times.magnitude, ep.times.magnitude+ep.durations.magnitude

    return starts, stops

def clean_trig(trig_times, t1,t2, anasig):
    cleaned = [ ]
    for i, trig_time in enumerate(trig_times):
        if is_artifact_in_window(anasig, trig_time+t1, trig_time+t2):
            print '   SKIP : artifact for trig', i, 'times ', trig_time
            continue
        cleaned.append(trig_time)
    return cleaned


def test_is_artifact_in_window():
    run = session.query(Run)[0]
    anasig = run.segments[0].analogsignals[7]

    print is_artifact_in_window(anasig,19.,20.)





def artifact_viewer():
    import pyqtgraph as pg


    app = pg.mkQApp()

    from OpenElectrophy.gui.viewers.multiviewer import MultiViewer
    view = MultiViewer()
    #~ run = session.query(Run)[0]
    #~ run = session.query(Run)[5]
    #~ run = Run.load(8, session = session)

    subject = session.query(Subject).filter_by(name = 'LEFC').first()
    run = subject.runs.filter_by(exp = 'R').filter_by(index = 3).first()


    neosigs = [ ] #a.to_neo() for a in run.segments[0].analogsignals ]


    epocharrays = [ ]
    for anasig in run.segments[0].analogsignals:
        name =  'Artifact {} {}'.format(anasig.name, anasig.id)
        ep = session.query(EpochArray).filter_by(name = name).first()
        if ep is None:
            ep = compute_and_store_artifact(anasig)
        neoep = ep.to_neo()

        color = mycolors.get_color(anasig)
        if color is not None:
            neosig = anasig.to_neo()
            neosig.annotate(color = color)
            neosigs.append(neosig)
            neoep.annotate(color = color)
            epocharrays.append(neoep)
    view.add_analogsignals( analogsignals = neosigs )
    view.add_epochs( epocharrays = epocharrays)

    view.timeseeker.set_start_stop(0, neosigs[0].t_stop.magnitude)

    for subviewer in view.subviewers:
        subviewer.viewer.xsize = 10.


    view.show()

    app.exec_()


if __name__ == '__main__':
    #test_detect_artifact_on_one_sig()
    test_compute_and_store_artifact()
    #~ test_is_artifact_in_window()
    #~ artifact_viewer()
