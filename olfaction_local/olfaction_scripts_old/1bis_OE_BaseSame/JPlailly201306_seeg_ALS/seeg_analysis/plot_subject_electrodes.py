# -*- coding: utf-8 -*-


import sys, os
sys.path.append('../behavior') 

from connection import *

from collections import OrderedDict

from PyQt4 import QtCore, QtGui
import pyqtgraph as pg

import pyplotbrain as ppb
import mycolors

import matplotlib.pyplot as plt

def hex_to_rgb(c1, alpha =1.):
    if c1 is None:
        return (.9, .9, .9, .9)
    c1 = str(c1)
    r = int(c1[1:3],16)/255.
    g = int(c1[3:5],16)/255.
    b = int(c1[5:7],16)/255.
    return (r, g, b, alpha)
    
def plot_one_subject(subject_name):
    subject = session.query(Subject).filter_by(name = subject_name).first()
    coords = [ ]
    colors = [ ]
    for rc in subject.blocks[0].recordingchannelgroups[0].recordingchannels:
        if rc.coordinate is None:
            #~ print rc.name, 'has no coordinate'
            continue
        coord = rc.coordinate.magnitude
        if np.any(np.isnan(coord)): continue
        
        color = mycolors.get_color(rc)
        if color is None:
            color = '#FFFFFF'
            print 'no color for', rc.name
        colors.append(color)
        coords.append(rc.coordinate)
    coords = np.array(coords)
    colors = np.array(colors)
    view = ppb.addView(with_config = True, cortical_alpha = .5)
    view.plot_mesh()
    
    for color in np.unique(colors):
        sel = colors==color
        if np.sum(sel)==0:continue
        view.add_node(coords[sel], color =  hex_to_rgb(color), size = 2)
    
    color =  (0,1,0,.8)
    view.add_node(np.array([[0,50,0]]), color = color, size = 4)    
    view.setWindowTitle(subject.name)
    return view


def plot_all_subject():
    views = [ ]
    for subject in session.query(Subject):
        views.append(plot_one_subject(subject.name))
    return views

if __name__ == '__main__':
    app = pg.mkQApp()
    
    #~ view = plot_one_subject('SEMC')
    #~ view = plot_one_subject('LEFC')
    #~ view = plot_one_subject('CHAF')
    #~ view = plot_one_subject('FERJ')
    #~ view = plot_one_subject('VACJ')
    
    views = plot_all_subject()
    
    app.exec_()
    
    