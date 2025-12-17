# -*- coding: utf-8 -*-

"""
Change this when OE3

"""

import quantities as pq
from datetime import datetime
import numpy as np


#~ import OpenElectrophy as OE
#~ print OE.classes
#~ from OpenElectrophy.core.base import OEBase
from OpenElectrophy.core import oeclasses, OEBase
from OpenElectrophy.gui.editdb import EditFieldsDialog

from collections import OrderedDict

import neo


kept = [
            'Block',
            'Segment',
            'AnalogSignal',
            'EpochArray',
            'EventArray',
            'RecordingChannel',
            'RecordingChannelGroup',
            ]

myclasses = [ ]
for c in oeclasses:
    if c.__name__ in kept:
        myclasses.append(c)
        
        # cleaan relationship
        for relname in ['one_to_many_relationship', 'many_to_one_relationship', 'many_to_many_relationship']:
            rel = getattr(c, relname)
            toremove = [ ]
            for e in rel:
                if e not in kept:
                    toremove.append(e)
            for e in toremove:
                rel.remove(e)
    #~ if c.__name__=='RecordingChannel':
        #~ c.attributes.append( ('color', str))
    




class Subject(OEBase):
    tablename = 'Subject'
    neoclass = None
    attributes =[ ('name', str),
                            ('num', str),
                            ('group', str),
                            ('encoding_order', str),
                            ('score_episodic_wwwS1', str),
                            ('score_episodic_wwwS4',str),
                            ('odor_pleasantness', np.ndarray),# pleasantness vector base 0
                            ]
    one_to_many_relationship = ['Run', 'Block', ]
    many_to_one_relationship = [  ]
    many_to_many_relationship = [ ]
    inheriting_quantities = None
    
    @property
    def recordingchannels(self):
        return self.blocks[0].recordingchannelgroups.filter_by(name = 'all channels')[0].recordingchannels


class Run(OEBase):
    tablename = 'Run'
    neoclass = None
    attributes =[  ('name', str),
                            ('index', int),
                            ('exp', str),
                            ('date', datetime),
                            ('group', str),
                            ('rand', str),
                            ('filename', str),
                            ('image', str),
                            ]
    one_to_many_relationship = [ 'Trial', 'Segment', 'EpochArray' ]#'RespirationSignal', 'EventArray', 'EpochArray',]
    many_to_one_relationship = [  ]
    many_to_many_relationship = [ ]
    inheriting_quantities = None


class Trial(OEBase):
    tablename = 'Trial'
    neoclass = None
    attributes =[         ('index', int),
                                    ('time', float),
                                    ('odor_num', int),  # base 1
                                    ('odor_name', str),
                                    ('wanted_odor_time', float),
                                    ('triggered_odor_time', float),
                                    ('first_inspiration_time', float),
                                    ('inspi_no_odor', float),
                                    ('recognition', int),
                                    ('recognition_time', float),
                                    ('context_time', float),
                                    ('selected_image', int),
                                    ('click_image_time', float),
                                    ('click_posistion', int),
                                    ('click_posistion_time', float),
                                    ('score_recognition', str),
                                    ('score_episodic_strict', str),
                                    ('score_episodic_elargi', str),
                                    ]
    one_to_many_relationship = [ ]
    many_to_one_relationship = [  ]
    many_to_many_relationship = [ ]
    inheriting_quantities = None



class RespirationSignal(OEBase):
    tablename = 'RespirationSignal'
    neoclass = None
    attributes =[('name', str),
                            ('channel_index', int),
                            ('t_start', pq.Quantity, 0),
                            ('sampling_rate', pq.Quantity, 0),
                            ('signal', pq.Quantity, 1),
                            ('cycle_times', pq.Quantity, 2),
                            ('odor_threshold', pq.Quantity, 0), 
                            ]
    one_to_many_relationship = [ ]
    many_to_one_relationship = [  'Segment',  'Run', ]
    many_to_many_relationship = [ ]
    inheriting_quantities = None


#~ class EventArray(OEBase):
    #~ tablename = 'EventArray'
    #~ neoclass = neo.EventArray
    #~ attributes =[ ('name', str),
                            #~ ( 'times', pq.Quantity, 1 ),
                            #~ ( 'labels',  np.ndarray,1 ,  np.dtype('S')) 
                                #~ ]
    #~ one_to_many_relationship = [ ]
    #~ many_to_one_relationship = [  ]
    #~ many_to_many_relationship = [ ]
    #~ inheriting_quantities = None

#~ class EpochArray(OEBase):
    #~ tablename = 'EpochArray'
    #~ neoclass = neo.EpochArray
    #~ attributes =[('name', str),
                            #~ ( 'times', pq.Quantity, 1 ),
                            #~ ( 'durations', pq.Quantity, 1 ),
                            #~ ( 'labels',  np.ndarray,1,  np.dtype('S'))  
                            #~ ]
    #~ one_to_many_relationship = [ ]
    #~ many_to_one_relationship = [  ]
    #~ many_to_many_relationship = [ ]
    #~ inheriting_quantities = None

#~ class AnalogSignal(OEBase):
    #~ tablename = 'AnalogSignal'
    #~ neoclass = neo.AnalogSignal
    #~ attributes =[('name', str),
                            #~ ( 'times', pq.Quantity, 1 ),
                            #~ ( 'durations', pq.Quantity, 1 ),
                            #~ ( 'labels',  np.ndarray,1,  np.dtype('S'))  
                            #~ ]
    #~ one_to_many_relationship = [ ]
    #~ many_to_one_relationship = [  ]
    #~ many_to_many_relationship = [ ]
    #~ inheriting_quantities = None


                                #~ ('signal', pq.Quantity, 1 ),
                                #~ ('sampling_rate', pq.Quantity, 0 ),
                                #~ ('t_start', pq.Quantity, 0 ),
                                #~ ('channel_index', int),
                                #~ ('name', str ),



#~ class EventArray(OEBase):
    #~ tablename = 'EventArray'
    #~ neoclass = None
    #~ attributes =[ ('name', str),
                            #~ ( 'times', pq.Quantity, 1 ),
                            #~ ( 'labels',  np.ndarray,1 ,  np.dtype('S')) 
                                #~ ]
    #~ one_to_many_relationship = [ ]
    #~ many_to_one_relationship = [  ]
    #~ many_to_many_relationship = [ ]
    #~ inheriting_quantities = None

#~ class EpochArray(OEBase):
    #~ tablename = 'EpochArray'
    #~ neoclass = None
    #~ attributes =[('name', str),
                            #~ ( 'times', pq.Quantity, 1 ),
                            #~ ( 'durations', pq.Quantity, 1 ),
                            #~ ( 'labels',  np.ndarray,1,  np.dtype('S'))  
                            #~ ]
    #~ one_to_many_relationship = [ ]
    #~ many_to_one_relationship = [  ]
    #~ many_to_many_relationship = [ ]
    #~ inheriting_quantities = None









#~ classes = [ Subject, Run, Trial, RespirationSignal, EventArray, EpochArray, AnalogSignal ]
classes = myclasses+ [ Subject, Run, Trial, RespirationSignal ]

for class_ in classes:
    class_.usable_attributes = OrderedDict()
    for attribute in class_.attributes:
        attrname, attrtype = attribute[0], attribute[1]
        class_.usable_attributes[attrname] = attrtype


class_by_name = { }
for class_ in classes:
    class_by_name[class_.__name__] = class_


# check bijectivity
for c1 in class_by_name.keys():
    for c2 in class_by_name[c1].one_to_many_relationship:
        if c1 not in class_by_name[c2].many_to_one_relationship :
            class_by_name[c2].many_to_one_relationship.append(c1)
    for c2 in class_by_name[c1].many_to_one_relationship:
        if c1 not in class_by_name[c2].one_to_many_relationship :
            class_by_name[c2].one_to_many_relationship.append(c1)
    for c2 in class_by_name[c1].many_to_many_relationship:
        if c1 not in class_by_name[c2].many_to_many_relationship :
            class_by_name[c2].many_to_many_relationship.append(c1)
    

