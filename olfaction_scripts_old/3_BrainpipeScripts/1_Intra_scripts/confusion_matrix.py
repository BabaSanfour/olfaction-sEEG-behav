import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from brainpipe.visual import *
import math

def plot_confusion_matrix(cm, xtickslabels1, xtickslabels2, ytickslabels, ylabel=None,xlabel=None,
                          normalize=False,cmap=plt.cm.Blues, size=(10,5), cbsides=1):
    """
    This function prints and plots the matrix
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')


    fig = plt.figure(figsize=size)
    ax1 = fig.add_subplot(111)
    ax1.set_axis_bgcolor("silver")

    if cbsides ==1:
        fmt = '.2f' #if normalize else 'd'
        vmin, vmax = np.nanmin(cm), np.nanmax(cm)
        thresh = np.nanpercentile(cm,80)
        for i, j in product(range(cm.shape[0]), range(cm.shape[1])):
            if not math.isnan(cm[i, j]):
                plt.text(j, i, format(cm[i, j], fmt),
                         horizontalalignment="center",verticalalignment="center",
                         color="black" if cm[i, j] < thresh else "white")
    else :
        #colorbar centered on zero
        cm_min, cm_max = np.nanmin(cm), np.nanmax(cm)
        print('min,max',cm_min,cm_max)
        vmin = cm_max*-1 if abs(cm_min) < cm_max else cm_min*-1
        vmax = cm_max if abs(cm_min) < cm_max else cm_min
        #color of data in the matrix
        fmt = '.2f' #if normalize else 'd'
        thresh1 = np.nanpercentile(cm,80)
        thresh2 = np.nanpercentile(cm,20)
        print(thresh1,thresh2)
        for i, j in product(range(cm.shape[0]), range(cm.shape[1])):
            if not math.isnan(cm[i, j]):
                plt.text(j, i, format(cm[i, j], fmt),
                         horizontalalignment="center",verticalalignment="center",
                         color="black" if thresh2 < cm[i, j] < thresh1 else "white")

    im = ax1.imshow(cm, interpolation='nearest', cmap=cmap, aspect='auto',vmin=vmin, vmax=vmax)
    fig.colorbar(im)
    nb = len(np.unique(xtickslabels1))
    #Deal with mutliple axes
    ytick_marks = np.arange(len(ytickslabels))
    xtick_marks1 = np.arange(nb/2,cm.shape[1], step=nb)
    xtick_marks2 = np.arange(cm.shape[1])
    
    ax1.set_xticks(xtick_marks1-0.5)
    ax1.set_xticklabels(xtickslabels2)
    ax1.set_yticks(ytick_marks)
    ax1.set_yticklabels(ytickslabels)
    
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(xtick_marks2)
    ax2.set_xticklabels(xtickslabels1)

    #Required to remove some white border
    ax1.autoscale(False)
    ax2.autoscale(False)
        
    #Add lines to separate freqs
    freq_marks = np.arange(cm.shape[1], step=nb)
    addLines(plt.gca(), vLines=freq_marks-0.5, vColor=['black']*len(freq_marks),
         vWidth=[2]*len(freq_marks), vShape=['-']*len(freq_marks))


    plt.tight_layout()
    plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    return fig