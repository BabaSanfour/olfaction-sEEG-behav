# -*- coding: utf-8 -*-


from connection import *

import pylab
import re
from scipy import stats

import neo
from OpenElectrophy.core import OEBase, neo_to_oe

import pandas

from align_bad_triggers import align_bad_triggers
from matplotlib import pyplot

subject_channel_trig = {
    'S00' : "G'8",   #'S00' : 'VII8',  #CHAF
    'S01' : 'MKR1+',  #VACJ
    'S02' : 'trig_Note', #SEMC
    'S03' : 'MKR1+',  #LEFC
    'S04' : 'trig_Note',  #FERJ
    'S05' : 'MKR1+',    #PIRJ
    'S06' : 'trig_Note',    #MICp
    }

channel_to_remove = {
    'CHAF' :  range(122,136),
    'VACJ' :  range(56,64) + range(124,130) + range(167,256),
    'SEMC' :  range(59,63) + range(123, 130),
    'LEFC' :  [63,128,129] + range(190,195)+range(219,256),
    'FERJ' :  range(127,130),
    'PIRJ' :  range(59,64) + range(126,136),
    'MICp' :  range(126,129), 
    }


verify_with_plot_bad_align = True


def signalToTrig(sig, thresh, front = True, clean = True, check = True):
    """
    Find pulse in a continuous signal.
    
    """
    if front:
        ind1, = np.where( (sig[:-1] <thresh) & (sig[1:]>=thresh) )
        ind2, = np.where( (sig[:-1] >=thresh) & (sig[1:]<thresh) )
    else:
        ind1, = np.where( (sig[:-1] >thresh) & (sig[1:]<=thresh) )
        ind2, = np.where( (sig[:-1] <=thresh) & (sig[1:]>thresh) )
    if clean:
        if ind1.size>=1 and ind2.size>=1 and ind2[0] < ind1[0]:
            ind2 = ind2[1:]
        if ind1.size>=1 and ind2.size>=1 and ind1[-1] > ind2[-1]:
            ind1 = ind1[:-1]
        if ind1.size == 1 and ind2.size==0:
            ind1= np.array([], dtype = int)
        if ind1.size == 0 and ind2.size==1:
            ind2= np.array([], dtype = int)
    if check:
        assert ind1.size == ind2.size, 'Not the same ind1 and ind2 size {} {}'.format(ind1.size, ind2.size)
    return ind1, ind2




def insert_one_run(filename_res, channel_trig_odor = None):

    filename_raw = filename_res.strip('.res')+'.raw'
    filename_header = filename_res.strip('.res')+'.header'
    filename_protocol = filename_res.strip('.res')+'.protocol'
    filename_param_olfacto = filename_res.strip('.res')+'.paraOlfacto'
    
    filename = os.path.basename(filename_res)
    assert filename.endswith('.res'), 'Filename do not ends with .res'
    filename = filename.strip('.res')
    
    subject_dir = os.path.dirname(os.path.dirname(filename_raw))
    #nomenclature
    info = { }
    for kv  in filename.split('_'):
        k,v = kv.split('=')
        info[k.lower()]= v
    
    subject_num = info['sujet']
    
    subject_name = os.path.basename(subject_dir).split('_')[1]

    
    
    dir_trc =  os.path.join(subject_dir, '{}_{}_SignalsSEEG'.format(subject_num, subject_name))

    filename_trc = os.path.join(dir_trc, filename+'.TRC')
    
    assert os.path.exists(filename_trc), 'Pas de fichier TRC {}'.format(filename_trc)
    
    
    
    
    subject = session.query(Subject).filter_by(name = subject_name).first()
    if subject is None:
        if info['exp'] == 'R':
            encoding_order = info['ordre']
            subject = Subject(name =subject_name,num = subject_num,  group = info['group'], encoding_order = info['ordre'])
        if info['exp'] == 'E':
            image = info['image']
            subject = Subject(name =subject_name,num = subject_num,  group = info['group'], image = info['image'])
        session.add(subject)
        session.commit()
    else:
        if info['exp'] == 'R':
            subject.encoding_order = info['ordre']
            session.commit()
        
        
    
    run = session.query(Run).filter_by(subject_id = subject.id).filter_by(index = info['run']).filter_by(exp = info['exp']).first()
    if run is not None:
        session.delete(run)
        session.commit()
    

    d ={ }
    for k in ['exp','group', 'rand', 'image']:
        if k in info:
            d[k] = info[k]
    d['index'] = info['run']
    l = tuple([ int(e) for e in info['date'].split('-') ])
    d['date'] = datetime.date(*l)
    
    run = Run(filename = filename, **d)
    subject.runs.append(run)
    session.commit()
    
    
    

        
    
    # read ephys files
    neo_bl = neo.MicromedIO(filename = filename_trc).read()[0]
    neo_seg = neo_bl.segments[0]
    neo_sigs = [ ]
    for neo_sig in neo_seg.analogsignals:
        if neo_sig.channel_index not in channel_to_remove[subject.name]:
            neo_sigs.append(neo_sig)
    neo_seg.analogsignals = neo_sigs
    
    if subject_name== 'CHAF':
        # patch pour les chiffres romain
        convertions = [ ('VIII' , 'H'), ('VII' , 'G'), ('VI' , 'F'), ('IV', 'D'), 
                                    ('IX' , 'I'), 
                                    ('IV', 'D'), ('III',  'C'), ('II' ,  'B'),
                                    ('X', 'J' ), ('V' ,  'E') , ('I' ,  'A'),
                                    ('8/', 'H' ),
                                    ]
        for neosig in neo_seg.analogsignals:
            #~ print neosig.name,
            for k, v in convertions:
                if k in neosig.name:
                    neosig.name = neosig.name.replace(k,v+"'")
                    break
            #~ print '>>>>',  neosig.name
    print
    if subject.blocks.count()==0:
        neo.io.tools.populate_RecordingChannel(neo_bl, remove_from_annotation = False)
        bl = neo_to_oe(neo_bl, dbinfo.mapped_classes, cascade = True)
        subject.blocks.append(bl)
        seg = bl.segments[0]
        for rc in bl.recordingchannelgroups[0].recordingchannels:
            rc.name = rc.analogsignals[0].name.lower().replace(' ' , '')
    else:
        seg = neo_to_oe(neo_seg, dbinfo.mapped_classes, cascade = True)
        rcg = subject.blocks[0].recordingchannelgroups.filter_by(name = 'all channels').first()
        for anasig in seg.analogsignals:
            rc = rcg.recordingchannels.filter_by(index = anasig.channel_index).first()
            rc.analogsignals.append(anasig)
    
    run.segments.append(seg)
    session.commit()

    #Read behavior files
    odors_num_to_name, odor_durations = open_param_olfacto(filename_param_olfacto)
    
    if run.exp == 'E':
        #results
        trials = open_results_encoding(filename_res)
    elif run.exp == 'R':
        #protocol and odors
        seq = open_protocol(filename_protocol)
        #results
        trials = open_results_recall(filename_res)
        
    # signals resp and trig odors on DACQ
    resp,  odors  = open_signals(filename_raw, filename_header)
    
    # regression between daq card and PC times : 
    txt_odor_times =  [ trial['triggered_odor_time'] for trial in trials ]
    a,b,r,tt,stderr=stats.linregress(txt_odor_times,odors.times.magnitude)
    print 'daq card and PC times  a', a, 'b', b, 'stderr', stderr
    for trial in trials:
        trial['time'] = trial['time']*a + b
        if 'wanted_odor_time' in trial:
            trial['wanted_odor_time'] = trial['wanted_odor_time']*a + b
        trial['triggered_odor_time'] = trial['triggered_odor_time']*a + b

    if  verify_with_plot_bad_align:
        fig, ax = pyplot.subplots()
        ax.set_title('txt file to odors (file RAW)')
        for i in txt_odor_times:
            ax.axvline(i, color = 'b')
        for i in odors.times.magnitude:
            ax.axhline(i, color = 'b')
        a,b,r,tt,stderr=stats.linregress(txt_odor_times,odors.times.magnitude)
        ax.plot(txt_odor_times, np.array(txt_odor_times)*a+b, color = 'r')
        ax.set_aspect('equal')



        

    
    run.respirationsignals.append(resp)
    run.epocharrays.append(odors)
    
    seg.respirationsignals.append(resp)
    seg.epocharrays.append(odors)
    
    seg.name = subject.name
    
    #~ session.commit()

    
    
    if channel_trig_odor is None:
        channel_trig = subject_channel_trig[subject.num]
        
    
    if  'trig' in channel_trig:
        _, name = channel_trig.split('_')
        ea = None
        for e in seg.eventarrays:
            if e.name == name:
                ea = e
                break
        micromed_tirg = neo.EventArray(name = 'Micromed_trig ({})'.format(channel_trig), times =ea.times)
        
    else:
        anasig_trig = seg.analogsignals.filter_by(name = channel_trig).first()
        assert anasig_trig is not None, 'No channel for trigger with {}'.format(channel_trig)
        sig = anasig_trig.signal.magnitude
        med = np.median(sig)
        max = np.max(sig)
        thresh = med+.7*(max-med)
        ind1, ind2 = signalToTrig(sig, thresh, front = True, clean = False, check = False)
        micromed_tirg = neo.EventArray(name = 'Micromed_trig ({})'.format(channel_trig), times = anasig_trig.to_neo().times[ind1])
    
    if subject.num== 'S01' and run.exp == 'E' and run.index == 2:
        # BUG ce jour de manip micromed a n'a pas le debut
        odors.times = odors.times[2:] 
        odors.durations = odors.durations[2:]
        trials = trials[2:]
        
        
    alignement_problem = False
    
    if micromed_tirg.times.size!=odors.times.size:
        print '   #####   @$&*%! WTF!  Not the same nb of trig :    {} vs {}'.format(odors.times.size, micromed_tirg.times.size)
        #~ print odors.times.magnitude
        #~ print micromed_tirg.times.magnitude
        #~ assert odors.times.size<micromed_tirg.times.size, 'more trigger on olfactometer than micromed ???? WTF really!'
        
        #~ print odors.times.magnitude
        #~ print micromed_tirg.times.magnitude
        mask = align_bad_triggers(odors.times.magnitude, micromed_tirg.times.magnitude,  thresh_stderr = 1e-3, thresh_slope = 1.e-2)
        #~ assert mask is not None, 'pas de solution pour les triggers'
        if mask is None:
            print '   #####   trigger NON resolved #####  WTF really!'
            alignement_problem = True
        else:
            print '   #####   trigger resolved merci sam (le vin bientot promis)'
        
        
        if  verify_with_plot_bad_align:
            fig, ax = pyplot.subplots()
            for i in odors.times.magnitude:
                ax.axvline(i, color = 'b')
            for i in micromed_tirg.times.magnitude:
                ax.axhline(i, color = 'c')
            for i in micromed_tirg.times.magnitude[mask]:
                ax.axhline(i, color = 'b')
            a,b,r,tt,stderr=stats.linregress(odors.times.magnitude, micromed_tirg.times.magnitude[mask])
            ax.plot(odors.times.magnitude, odors.times.magnitude*a+b, color = 'r')
            ax.set_aspect('equal')
            ax.set_title(filename_res)
        
        micromed_tirg.times = micromed_tirg.times[mask]
    
    if  verify_with_plot_bad_align:
        fig, ax = pyplot.subplots()
        ax.set_title('odors (file RAW) to micromed')
        for i in odors.times.magnitude:
            ax.axvline(i, color = 'b')
        for i in micromed_tirg.times.magnitude:
            ax.axhline(i, color = 'b')
        a,b,r,tt,stderr=stats.linregress(odors.times.magnitude, micromed_tirg.times.magnitude)
        ax.plot(odors.times.magnitude, odors.times.magnitude*a+b, color = 'r')
        ax.set_aspect('equal')
        
    
    
    #~ assert micromed_tirg.times.size==odors.times.size, 'Not the same nb of trig {} {}'.format(micromed_tirg.times.size, odors.times.size)
    seg.eventarrays.append(neo_to_oe(micromed_tirg,dbinfo.mapped_classes, cascade = True))
    

    if not alignement_problem:
        # times corrections of daq card triggers with micromed files:
        a,b,r,tt,stderr=stats.linregress(odors.times.magnitude, micromed_tirg.times.magnitude)
        
        if (subject.num== 'S00' and run.exp == 'E' and run.index == 1) or  \
             (subject.num== 'S05' and run.exp == 'E' and run.index == 1):
            #Le seigneur nous a envoyé une epreuve : le trig numero 11/5 est décalé.
            # Personne ne sait pourquoi, donc onle virede l aregression uniqement.
            to_remove = {'S00' : 11, 'S05': 3}[subject.num]
            vect1, vect2 = odors.times.magnitude.copy(), micromed_tirg.times.magnitude.copy()
            vect1[to_remove] = np.nan
            vect2[to_remove] = np.nan
            a,b,r,tt,stderr=stats.linregress(vect1[~np.isnan(vect1)], vect2[~np.isnan(vect2)])
            
        
        print 'daq card on TRC a', a, 'b', b, 'stderr', stderr
        odors.times = odors.times* a + b*pq.s
        resp.sampling_rate /= a
        resp.t_start += b*pq.s
    
        # times corrections  PC times with micromed files: 
        for trial in trials:
            trial['time'] = trial['time']*a+b
            trial['triggered_odor_time'] = trial['triggered_odor_time']*a+b
            if 'wanted_odor_time' in trial:
                trial['wanted_odor_time'] = trial['wanted_odor_time']*a+b
        
    
    # complete info
    for t, trial in enumerate(trials):
        trial['index'] = t
        #~ print trial['odor_num']
        #~ print odors_num_to_name
        trial['odor_name'] =  odors_num_to_name[trial['odor_num']]
        trial['first_inspiration_time'] = None
    
    
    # create
    for t, trial in enumerate(trials):
        t = Trial(**trial)
        run.trials.append(t)
    session.commit()
    
    if alignement_problem:
        #Remove signal because they are not align with respiration
        run.segments.remove(seg)
        session.commit()
    
    assert not alignement_problem, 'problem de nombre de trig'
    
    
        
        

def test_insert_one_run():
    #~ dir = '/home/samuel/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir = 'N:/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'

    
    #~ dir +='S00_CHAF_ordre=FM/S00_CHAF_Behavior/'
    #~ filename_res = dir+'Date=2013-06-25_Sujet=S00_Exp=E_Run=1_Group=1_Image=F.res'
    
    
    #~ dir +='S02_SEMC_ordre=FM/S02_SEMC_Behavior/'
    #~ filename_res = dir+'Date=2013-10-25_Sujet=S02_Exp=R_Run=1_Group=1_Rand=1_Ordre=FM.res'
    #~ filename_res = dir+'Date=2013-10-24_Sujet=S02_Exp=E_Run=2_Group=1_Image=M.res'
    
    
    #~ dir +='S01_VACJ_ordre=MF/S01_VACJ_Behavior/'
    #~ filename_res = dir+'Date=2013-10-18_Sujet=S01_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.res'
    #~ filename_res = dir+'Date=2013-10-17_Sujet=S01_Exp=E_Run=2_Group=1_Image=F.res'
    
    #~ dir +='S02_SC_ordre=FM/S02_SC_Behavior/'
    #~ filename_res = dir+'Date=2013-10-25_Sujet=S02_Exp=R_Run=2_Group=1_Rand=1_Ordre=FM.res'
    #~ filename_res = dir+'Date=2013-10-23_Sujet=S02_Exp=E_Run=1_Group=1_Image=F.res'

    #~ dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    #~ filename_res = dir+'Date=2013-11-29_Sujet=S03_Exp=R_Run=2_Group=1_Rand=1_Ordre=MF.res'


    #~ dir +='S05_PIRJ_ordre=MF/S05_PIRJ_Behavior/'
    #~ filename_res = dir+'Date=2014-12-03_Sujet=S05_Exp=E_Run=1_Group=1_Image=M.res'
    #~ filename_res = dir+'Date=2014-12-04_Sujet=S05_Exp=E_Run=2_Group=1_Image=F.res'

    dir +='S06_MICp_ordre=FM/S06_MICp_Behavior/'
    #~ filename_res = dir+'Date=2015-10-12_Sujet=S06_Exp=E_Run=1_Group=1_Image=F.res'
    #~ filename_res = dir+'Date=2015-10-13_Sujet=S06_Exp=E_Run=2_Group=1_Image=M.res'
    #~ filename_res = dir+'Date=2015-10-14_Sujet=S06_Exp=R_Run=1_Group=1_Rand=1_Ordre=FM.res'
    #~ filename_res = dir+'Date=2015-10-14_Sujet=S06_Exp=R_Run=2_Group=1_Rand=1_Ordre=FM.res'
    #~ filename_res = dir+'Date=2015-10-14_Sujet=S06_Exp=R_Run=3_Group=1_Rand=1_Ordre=FM.res'

    
    
    insert_one_run(filename_res)

    

    

def insert_all_run():
    #~ dirname = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    #~ dirname = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/S00_CHAF_ordre=FM/'
    #~ dirname = 'D:/EPISODIC/Functional/seeg_episodic\DATA_BEHAVIOR_FUNCTIONAL'
    
    #~ dirname = '/media/data/19_EPISODIC_SEEG/DATA_BEHAVIOR_FUNCTIONAL/'
    #~ dirname = '/home/alsaive/smb4k/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    #~ dirname = '/home/samuel/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/S06_MICp_ordre=FM/'
    #~ dirname = 'N:/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/S06_MICp_ordre=FM/'
    #~ dirname = 'S06_MICp_ordre=FM/'
    #~ dirname = 'S02_SEMC_ordre=FM/'
    
    #la totale
    #~ dirname = '/home/samuel/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL'
    dirname = '/home/karim/Documents/Dropbox/Intra_EM/0_RawData_trc_avi_behav/S06_MICp_ordre=FM/'
    
    errors = [ ]
    for root, dirs, files in os.walk(dirname):
        for file in files:
            if file.endswith('.res'):
                print file
                
                try:
                #~ if True:
                    filename_res = os.path.join(root, file)
                    insert_one_run(filename_res)
                except AssertionError as e:
                    print '    !! ERREUR ', e

                    errors.append((os.path.basename(file),e))
                except :
                    print '    !! ERREUR inconnue'
                    errors.append((os.path.basename(file), 'inconnue'))
                #~ print

    print '#### ERRORS ####'
    for e,t in errors:
        print e, t
    print ' nb errors = ', len(errors)



def open_results_recall(filename_res_recall):
    
    trials = [ ]
    for line in open(filename_res_recall, 'rU'):
        #~ print line
        # TODO
        v = line.strip('\n').split('\t')
        #~ print 'v[4]', v[4]
        #~ print v
        trial = dict(time = float(v[0])/1000.,
                        odor_num = int(v[1]),
                        wanted_odor_time = float(v[2])/1000.,
                        triggered_odor_time = float(v[3])/1000.,
                        recognition = {'1':1, '2' : 0, '0':None,}[v[4]],
                        recognition_time = float(v[5])/1000.,
                        context_time = float(v[6])/1000.,
                        click_posistion = int(v[7]),
                        click_posistion_time = float(v[8])/1000.,
                        selected_image = int(v[9]),
                        click_image_time = float(v[10])/1000.,
                        )
        #~ print trial
        trials.append( trial )
    return trials
    
def test_open_results_recall():
    dir = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    filename_res_recall = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.res'

    trials = open_results_recall(filename_res_recall)
    for t in trials:
        print t

def open_results_encoding(filename_res_encoding):
    
    trials = [ ]
    for line in open(filename_res_encoding, 'rU'):
        #~ print line
        # TODO
        v = line.strip('\n').split('\t')
        #~ print 'v[4]', v[4]
        #~ print v
        trial = dict(
                        odor_num = int(v[0]),
                        time = float(v[1])/1000.,
                       triggered_odor_time = float(v[2])/1000.,
                        )
        #~ print trial
        trials.append( trial )
    return trials
    
def test_open_results_encoding():
    dir = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    filename_res_encoding = dir+'Date=2013-11-27_Sujet=S03_Exp=E_Run=1_Group=1_Image=M.res'

    trials = open_results_encoding(filename_res_encoding)
    for t in trials:
        print t


def open_protocol(filename_protocol):
    seq = [ ]
    for line in open(filename_protocol, 'rU'):
        #~ print line
        seq.append(int(line.strip('\n')))
    return seq
    
def test_open_protocol():
    dir = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    filename_protocol = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.protocol'
    
    seq = open_protocol(filename_protocol)
    print seq


def open_param_olfacto(filename_param_olfacto):
    odors_num_to_name = { }
    odor_durations = { }
    for line in open(filename_param_olfacto, 'rU'):
        line = line.decode('iso-8859-15')
        r = re.findall('(\d+) (\d+) (\d+) (\d+) (\d+) ([\S ]+)', line)
        if len(r) == 1: 
            num, _,_,_, odor_duration, odor_name = r[0]
            #~ print line.strip('\n').split(' ')
            #~ num, _,_,_, odor_duration, odor_name = line.strip('\n').split('\t')
            num = int(num)
            odors_num_to_name[num] = odor_name 
            odor_durations[num] = float(odor_duration)/1000.*pq.s
    
    return odors_num_to_name, odor_durations



    
def test_open_param_olfacto():
    dir = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    filename_header = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.header'
    filename_param_olfacto = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.paraOlfacto'


    odors_num_to_name, odor_durations = open_param_olfacto(filename_param_olfacto)
    print odors_num_to_name
    print odor_durations








def open_signals(filename_raw, filename_header):
    
    f = open(filename_header, 'rU')
    sampling_rate = float(f.readline().strip('frequence:').strip('\n'))*pq.Hz
    dt = np.dtype(f.readline().strip('dtype:').strip('\n'))
    nb = int(f.readline().strip('nbvoies:').strip('\n'))
    channels = [ ]
    for i in range(nb):
        channels.append(f.readline().strip('nom%d:'%(i+1)).strip('\n'))
    
    sigs = fromfile(filename_raw, dtype = dt).reshape(-1, len(channels))
    t_vect = (arange(0, sigs.shape[0]).astype('f')/sampling_rate).rescale('s')
    #~ print t_vect

    
    resp, odor_threshold, odors, trig_irm = (None,)*4
    t_start = 0*pq.s
    for i, channel in enumerate(channels):
        sig=  sigs[:,i]
        if channel == 'Respi':
            resp = RespirationSignal( signal = sig*pq.V, name = 'respiration', t_start = 0.*pq.s, 
                                    sampling_rate = sampling_rate)
        
        if channel == 'Seuil':
            odor_threshold = mean(sig)*pq.V
            
        if channel == 'Odeur':
            thresh = sig.min()+ (sig.max()-sig.min())*.7
            ind1, ind2 = signalToTrig(sig, thresh, front = True, clean = True, check = True)
            starts = t_vect[ind1]
            stops = t_vect[ind2]
            odors = EpochArray(times = starts, durations = stops-starts , name = 'odors')

    #~ print 't_start', t_start
    resp.odor_threshold = odor_threshold
    resp.t_start = t_start
    odors.times = odors.times + t_start
    
    return resp,  odors


    
def test_open_signals():
    dir = '/home/sgarcia/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/DATA_BEHAVIOR_FUNCTIONAL/'
    dir +='S03_LEFC_ordre=MF/S03_LEFC_Behavior/'
    filename_header = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.header'
    filename_raw = dir+'Date=2013-11-29_Sujet=S03_LEFC_Exp=R_Run=1_Group=1_Rand=1_Ordre=MF.raw'
    
    resp,  odors = open_signals(filename_raw, filename_header)
    
    times = resp.t_start + np.arange(resp.signal.shape[0])/resp.sampling_rate
    print times.shape
    
    
    print odors.times.size
    
    pylab.plot(times, resp.signal)
    pylab.axhline(resp.odor_threshold, color ='g')
    for t in odors.times:
        pylab.axvline(t, color ='r')
    
    pylab.show()



#~ def insert_hedonicity():
    #~ filename = 'pleasantness_raw_data.csv'
    #~ df = pandas.read_csv(filename, sep = '\t', header = 0, na_values = 'None')
    #~ print df
    
    #~ for subject in session.query(Subject):
        #~ if subject.name in df.columns:
            #~ subject.odor_pleasantness = df[subject.name].view(np.ndarray )/1.5
            #~ subject.save(session = session)
        #~ else:
            #~ print 'pas', subject.name
    
    

def insert_one_electrodes_coodinate(filename):
    _, _, subject_name, _,_ = os.path.basename(filename).split('_')
    
    subject_name = subject_name.upper()
    print subject_name   
    
    subject = session.query(Subject).filter_by(name = subject_name).first()
    if subject is None:
        print subject_name, ' is not in DB'
        return
    
    filename2 = filename.replace('Pos.txt', 'Name.txt')
    names = np.loadtxt(filename2, delimiter = '\t', dtype = str)
    names = np.array([ name.replace('p', "'").lower()  for name in names ], dtype = str)

    
    try:
        coords = np.loadtxt(filename, delimiter = '\t')
    except:
        coords = np.loadtxt(filename, delimiter = ' ')
    
    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels:
        rc_name = rc.name.replace(' ', '').lower()
        print rc_name
        ind,  = np.where(names==rc_name)
        if len(ind)==1:
            rc.coordinate = coords[ind[0], :] * pq.cm
            session.commit()
        else:
            print 'not in names', rc_name
    

def insert_all_electrodes_coordinate():
    
    #~ dirname = 'S02_SEMC_ordre=FM/'

    #~ dirname = '/home/samuel/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/Infos_Patients'
    dirname = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI'
    #~ dirname = 'D:/EPISODIC/Functional/seeg_episodic/Infos_Patients'
    
    errors = [ ]
    for root, dirs, files in os.walk(dirname):
        if 'old' in root.lower():
            continue
        for file in files:
            if not file.endswith('MNI_Pos.txt'): continue
            print file
            try:    
                insert_one_electrodes_coodinate(os.path.join(root, file))
            except:
                print 'erreur'


def create_channel_groups():
    
    for bl in session.query(Block):
        for rcg in bl.recordingchannelgroups:
            if rcg.name  != 'all channels':
                session.delete(rcg)
                session.commit()
        rcg0 =  bl.recordingchannelgroups[0]
        short_names = [ ]
        for rc in rcg0.recordingchannels:
            #~ print rc.name
            short_name = str(rc.name)
            for l in '0123456789':
                short_name = short_name.replace(l, '')
            short_names.append(short_name)
        short_names = np.array(short_names)
        
        for group_name in np.unique(short_names):
            if '+' in group_name: continue
            if '@' in group_name: continue
            if '/' in group_name: continue
            if '.' in group_name: continue
            sel = group_name==short_names
            if np.sum(sel)<=1: continue
            rcg = RecordingChannelGroup(name = 'Group {}'.format(group_name))
            bl.recordingchannelgroups.append(rcg)
            for ind in np.where(sel)[0]:
                rcg.recordingchannels.append(rcg0.recordingchannels[ind])
            session.commit()
        for rcg in bl.recordingchannelgroups:
            print rcg.name
        print 


def insert_electrode_labels():
    
    #~ dirname = '/home/samuel/mnt/CRNLDATA/crnldata/projets_communs/seeg_episodic/Infos_Patients/'
    dirname = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/'
    
    for root, dirs, files in os.walk(dirname):
        for file in files:
            if file.endswith('.xlsx') and 'Labellisation' in file:
                print file
                filename = os.path.join(root, file)
                try:
                #if True:
                    #~ filename = '/home/sgarcia/LabellisationCHAF.xlsx'
                    df = pandas.io.excel.read_excel(filename,0)
                    #print (df['label'])
                    #df['plot'] = df['plot'].map( lambda s: s.lower())
                    #df['plot'] = df['plot'].map( lambda s: s.lower())

                    subject_name = os.path.basename(filename).replace('Labellisation','').replace('.xlsx','')
                    print subject_name
                    subject = session.query(Subject).filter_by(name = subject_name).first() 
                    #~ print subject
                    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels:
                        #~ print rc
                        sel = rc.name == df['plot']
                        if np.sum(sel)==1:
                            rc.description =df['label'][sel].tolist()[0]
                        #~ print rc
                    session.commit()
                except:
                    print 'ça merde'
    
    

    #~ print df.columns
    #~ print df

if __name__ =='__main__':
    """
    Pour creer base vierge:
      * insert_all_run()
      * insert_all_electrodes_coordinate()
      * create_channel_groups()
      * lancer score_encoding, score _reconnaissance, score_episodic
      * insert_electrode_labels
    
    """
    #~ test_open_signals()
    #~ test_open_protocol()
    #~ test_open_results_recall()
    #~ test_open_results_encoding()
    
    #~ test_open_param_olfacto()
    
    #test_insert_one_run()
    
    #insert_all_run()
    
    #~ insert_hedonicity()
    
    #insert_all_electrodes_coordinate()
    
    create_channel_groups()
    
    
    #insert_electrode_labels()
    
    #if verify_with_plot_bad_align:
    #    pyplot.show()




