# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 
from connection import *

#~ from sqlalchemy import inspect


electrodegroup_dict = {
        'a' : ('Amygdala', '#FFBF00'),
        'b':  ('Hippocampus', '#FF0000'),
        'd':  ('Entorhinal area', '#FF8000'),
        'e':  ('Fronto-polar/Temporal2', '#31B404'),
        'f':  ('Fronto-polar', '#01DF3A'),
        'g':  ('post Cingulate','#2EC5FF'),
        'h':  ('Temporal1 post',  '#F7D358'),
        'j':  ('Temporal pole', '#F5BCA9'),
        'k':  ('ant Cingulate', '#FFFF00'),
        'l':  ('Temporal3 post', '#FA5858'),
        'm':  ('Precentral', '#81F79F'),
        'o':  ('Orbitofrontal', '#8000FF'),
        'q':  ('Parietal sup', '#2E64FE'),
        'r':  ('Parietal inf',  '#0101DF'),
        's':  ('Supl motor area',  '#DF013A'),
        't':  ('Temporal1 middle', '#FE642E'),
        'u':  ('Frontal3 middle', '#86B404'),
        'v':  ('Occipital', '#DF01D7'),
        'w':  ('Parietal/occipital', '#084B8A'),
        'x':  ('Parietal', '#7401DF'),
        'y':  ('Angular', '#B404AE'),
        'z':   ('Frontal1_2 middle',  '#04B486'),
        }



def get_color_RecordingChannel(rc):
    color = None
    if  rc.name[1] in "'0123456789":
        l = rc.name[0]
        description, color = electrodegroup_dict.get(l, (None, None))
    return color

def get_color_AnalogSignal(anasig):
    rc = anasig.recordingchannel
    return get_color_RecordingChannel(rc)


def get_color(obj):
    if obj.__class__.__name__=='AnalogSignal':
        return get_color_AnalogSignal(obj)
    if obj.__class__.__name__=='RecordingChannel':
        return get_color_RecordingChannel(obj)
        



def reset_all_color():
    for rc in session.query(RecordingChannel):
        rc.color = None
        rc.description = None
        session.commit()




if __name__ == '__main__':
    reset_all_color()