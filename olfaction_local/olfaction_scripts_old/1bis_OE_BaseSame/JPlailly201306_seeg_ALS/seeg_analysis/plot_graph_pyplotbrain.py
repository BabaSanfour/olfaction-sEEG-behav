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
    
    
def return_electrode_colors(channel_names,channel_coords):
    
    display_channel = np.zeros(shape = (len(channel_names)),dtype = 'bool')
    
    display_channel_colors = []
    
    for i,name in enumerate(channel_names):
        
        #print name
        #print channel_coords[i,:]
        
        if np.any(np.isnan(channel_coords[i,:])): continue
        
        l = name[0]
        
        description, color = mycolors.electrodegroup_dict.get(l, (None, None))
        
        if color is None:
            color = '#FFFFFF'
            print 'no color for', name
            
        #print color
        
        display_channel_colors.append(color)
        
        display_channel[i] = True
    
    print np.sum(display_channel == True)/float(display_channel.shape[0])
    
    colors = np.array(display_channel_colors)
    
    #print colors
    
    return display_channel,colors
        
def plot_electrodes(channel_names,channel_coords):
    
    display_channel,colors = return_electrode_colors(channel_names,channel_coords)
    
    coords = channel_coords[display_channel,:]
    
    print coords
    
    #TODO FIXME: pourquoi Y est a l'envers ?????
    #~ coords[:,1] = -coords[:,1]
    
    
    app = pg.mkQApp()
    view = ppb.addView(with_config = True)#, cortical_alpha = .4)
    
    view.plot_mesh()
    
    
    for color in np.unique(colors):
        sel = colors==color
        if np.sum(sel)==0:continue
        #~ print color, np.sum(sel)
        view.add_node(coords[sel], color =  hex_to_rgb(color), size = 2)
    
    color =  (0,1,0,.8)
    view.add_node(np.array([[0,50,0]]), color = color, size = 4)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    
    import matplotlib.colors as colors

    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),cmap(np.linspace(minval, maxval, n)))
    
    return new_cmap

def plot_electrodes_graphs(channel_names,channel_coords,signif_diff_mat,export_path,pref):
    
    display_channel,colors = return_electrode_colors(channel_names,channel_coords)
    
    coords = channel_coords[display_channel,:]
    
    #TODO FIXME: pourquoi Y est a l'envers ?????
    #~ coords[:,1] = -coords[:,1]
    
    
    app = pg.mkQApp()
    view = ppb.addView(with_config = True)#, cortical_alpha = .4)
    
    view.plot_mesh()
    
    
    for color in np.unique(colors):
        sel = colors==color
        if np.sum(sel)==0:continue
        #~ print color, np.sum(sel)
        view.add_node(coords[sel], color =  hex_to_rgb(color), size = 2)
    
    #color =  (0,1,0,.8)
    #view.add_node(np.array([[0,50,0]]), color = color, size = 4)
    
    ################## adding connectivity matrix
    
    scale_signif_diff_mat = signif_diff_mat[display_channel,:][:,display_channel]
    
    print scale_signif_diff_mat.shape
    
    cmap = plt.get_cmap('jet')
    
    cmap_vals = cmap(np.linspace(0.2,0.8,9))
    
    print cmap_vals
    
    #0/0
    
    
    #cmap_vals = truncate_colormap(cmap, 0.1, 1.0)(np.arange(9))
    
    ##cmap_vals = np.array(plt.get_cmap('spectral',9)(np.arange(9)))
    

    
    
    
    #print cmap_vals
    
    #print np.where(scale_signif_diff_mat == -4)
    ##0/0
    
    #cmap_vals[0,:] = (.1,.1,.1,1.0)
    
    cmap_vals[:,3] = .7
    
    
    view.add_edge(coords,np.array(scale_signif_diff_mat == 1,dtype = int) * 2,color = cmap_vals[5])
    view.add_edge(coords,np.array(scale_signif_diff_mat == 2,dtype = int) * 4,color = cmap_vals[6])
    view.add_edge(coords,np.array(scale_signif_diff_mat == 3,dtype = int) * 8,color = cmap_vals[7])
    view.add_edge(coords,np.array(scale_signif_diff_mat == 4,dtype = int) * 16,color = cmap_vals[8])
    
    view.add_edge(coords,np.array(scale_signif_diff_mat == -1,dtype = int) * 2,color = cmap_vals[3])
    view.add_edge(coords,np.array(scale_signif_diff_mat == -2,dtype = int) * 4,color = cmap_vals[2])
    view.add_edge(coords,np.array(scale_signif_diff_mat == -3,dtype = int) * 8,color = cmap_vals[1])
    view.add_edge(coords,np.array(scale_signif_diff_mat == -4,dtype = int) * 16,color = cmap_vals[0])
    
    view.resize(1200,800)
    
    ## from right
    view.glview.setCameraPosition(distance = 250,elevation = 0, azimuth = 0)
    view.to_file(os.path.join(export_path,pref + '_from_right.png'))
    
    ### from front
    view.glview.setCameraPosition(distance = 280,azimuth = 90,elevation = 0)
    view.to_file(os.path.join(export_path,pref + '_from_front.png'))

    ### from top
    view.glview.setCameraPosition(distance = 250,azimuth = 180,elevation = 90)
    view.to_file(os.path.join(export_path,pref + '_from_top.png'))

    ### from left
    view.glview.setCameraPosition(distance = 250,azimuth = 180,elevation = 0)
    view.to_file(os.path.join(export_path,pref + '_from_left.png'))

    #self.glview .setCameraPosition(160,160,15)
    app.exec_()
    
if __name__ == '__main__':
    
    test1()
    