## \package MD.analysis.pyplot 
# \brief This module is used for making simple plots on the fly using matplotlib

import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def params_small(width=117,mult=1,xmult=1):
    from math import sqrt
    #Set Latex Parameters 
    fig_width_pt = width*xmult  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean*mult      # height in inches
    fig_size =  [fig_width,fig_height]
    mpl.rc('font', family='helvet')
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'text.fontsize': 8,
              'legend.fontsize': 5,
              'xtick.labelsize': 7,
              'ytick.labelsize': 7,
              'text.usetex': True,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
def params(width=235,mult=1):
    from math import sqrt
    #Set Latex Parameters 
    fig_width_pt = width  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height*mult]
    mpl.rc('font', family='serif')
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 6,
              'xtick.labelsize': 7,
              'ytick.labelsize': 7,
              'text.usetex': True,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
#save current figure in certaint format
def save_fig(save='fig',savefmt='png',fig_size=[3.251,2.00]):
    plt.savefig(save+'.'+savefmt,figsize=fig_size,dpi=300,
            bbox_inches='tight',format=savefmt,orientation='portrait',
            pad_inches=0.03)
#### The Main Plotting function
def plot(X, Y, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, y, marker=marker, linestyle=linestyle, label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if showleg:
        ax1.legend(loc=2)
    save_fig(save=save,fig_size=fig_size)
def plot2(X1,Y1,X2,Y2,xlabel='x',ylabel='y',label1='',
        label2="",save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1,Y1,'x-',label=label1)
    ax1.plot(X2,Y2,'s-',label=label2)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if showleg:
        ax1.legend(loc=2)
    save_fig(save=save,fig_size=fig_size)
def plot3(X,Y1,Y2,Y3,xlabel='x',ylabel='y',label1='',
        label2="",label3="",save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X,Y1,color = '#401624',label=label1)
    ax1.plot(X,Y2,color = '#254117',label=label2)
    ax1.plot(X,Y3,color = '#162440',label=label2)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if showleg:
        ax1.legend(loc=1)
    save_fig(save=save,fig_size=fig_size)
def plot_multi(X, Y, xlabel='x', ylabel='y',save='fig'):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for i in range(y.shape[0]):
        ax1.plot(x, y[i])
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    save_fig(save=save,fig_size=fig_size)
#### histogram plottiong
def plot_bar(X,Y,xlabel='x',ylabel='y',title='',label='',
        save='fig1',showleg=False,savefmt='png'):
    #Plot Data
    fig_size=params()
    y_array = np.array(Y)
    x_array = np.array(X)
    fig1=plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.bar(x_array,y_array,0.0025)
    if showleg:
        ax1.legend(loc=2)
    ax1.set_xlabel(xlabel,labelpad=2)
    ax1.set_ylabel(ylabel,labelpad=2)
    ax1.axes.xaxis.get_major_locator()._nbins=7 
    ax1.axes.yaxis.get_major_locator()._nbins=5 
    save_fig(fig_size=fig_size,save=save,savefmt='png')

