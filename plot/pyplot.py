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
              'legend.fontsize': 4,
              'xtick.labelsize': 7,
              'ytick.labelsize': 7,
              'text.usetex': False,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
def params_square(width=235,mult=1):
    from math import sqrt
    #Set Latex Parameters 
    fig_width_pt = width  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width                  # height in inches
    fig_size =  [fig_width,fig_height]
    mpl.rc('font', family='serif')
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 4,
              'xtick.labelsize': 7,
              'ytick.labelsize': 7,
              'text.usetex': False,
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
        label='', save='fig', showleg=False,
        limx=False,limy=False,logx=False,logy=False):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params_square()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, y, linestyle, marker=marker, linestyle=linestyle, label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if logx:
        ax1.set_xscale('log')
    if logy:
        ax1.set_yscale('log')
    if limx != False:
        plt.xlim(limx)
    if limy != False:
        plt.ylim(limy)
    if showleg:
        ax1.legend(loc=2)
    if save != False:
        save_fig(save=save,fig_size=fig_size)
def plot2(X1,Y1,X2,Y2,xlabel='x',ylabel='y',label1='',
        label2='',save='fig',showleg=False,limx=False,limy=False,logx=False,logy=False, 
        show = False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1,ls='',marker='x', color = '#990000', label=label1)
    ax1.plot(X2, Y2, color = '#336622', label=label2)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if logx:
        ax1.set_xscale('log')
    if logy:
        ax1.set_yscale('log')
    if limx != False:
        plt.xlim(limx)
    if limy != False:
        plt.ylim(limy)
    if showleg:
        ax1.legend(loc=2)
    if save != False:
        save_fig(save=save,fig_size=fig_size)
    if show:
        plt.show()
    save_fig(save=save,fig_size=fig_size)
def plot3(X1,Y1,X2,Y2,X3,Y3,xlabel='x',ylabel='y',label1='',
        label2="",label3="",save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000', label=label1)
    ax1.plot(X2, Y2, color = '#CC9900', label=label2)
    ax1.plot(X3, Y3, color = '#0099FF', label=label3)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if lim != False:
        plt.xlim(lim)
    if showleg:
        ax1.legend(loc=1)
    fig1.show()
    save_fig(save=save,fig_size=fig_size)
def plot4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4='',save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000', label=label1)
    ax1.plot(X2, Y2, color = '#CC9900', label=label2)
    ax1.plot(X3, Y3, color = '#0099FF', label=label3)
    ax1.plot(X4, Y4, color = '#006666', label=label4)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #plt.xlim([0,2.2e8])
    #plt.ylim([0,1])
    if showleg:
        ax1.legend(loc=4)
    save_fig(save=save,fig_size=fig_size)
def plot_multi(X, Y, label, xlabel='x', ylabel='y',
        save='fig',title=False,split=False,limx=False,
        limy=False,logx=False,logy=False,loc=3,showleg=False,
        make_marker=False,M=['x','>','s','o','d']):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params_square()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for i in range(y.shape[0]):
        if split==True:
            if i%2 == 1:
                ax1.plot(x[i], y[i], label = label[i])
            else:
                ax1.plot(x[i], y[i], label = label[i],linestyle='',marker='x')
        elif make_marker:
            ax1.plot(x[i], y[i], label = label[i],linestyle='',marker=M[i])
        else:
            ax1.plot(x[i], y[i], label = label[i])
    ax1.legend(loc=loc)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if title != False:
        plt.title(title)
    if logx:
        ax1.set_xscale('log')
    if logy:
        ax1.set_yscale('log')
    if limx != False:
        plt.xlim(limx)
    if limy != False:
        plt.ylim(limy)
    if showleg:
        ax1.legend(loc=loc)
    if save != False:
        save_fig(save=save,fig_size=fig_size)
#### histogram plottiong
def plot_bar(X,Y,xlabel='x',ylabel='y',title='',label='',
        save='fig1', showleg=False, savefmt='png',
        width=0.1,color='b',xmax=4,ymax=1):
    #Plot Data
    fig_size=params()
    y_array = np.array(Y)
    x_array = np.array(X)
    fig1=plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.bar(x_array,y_array,width,color=color)
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    if showleg:
        ax1.legend(loc=2)
    ax1.set_xlabel(xlabel,labelpad=2)
    ax1.set_ylabel(ylabel,labelpad=2)
    ax1.axes.xaxis.get_major_locator()._nbins=7
    ax1.axes.yaxis.get_major_locator()._nbins=5
    if save != False:
        save_fig(fig_size=fig_size,save=save,savefmt='png')
def plot_msd_int(X, Y, linestyle='-', marker='', xlabel='x', ylabel='y',
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
    plt.xlim([1.3e8,2.4e8])
    if showleg:
        ax1.legend(loc=2)
    save_fig(save=save,fig_size=fig_size)
def plot_sf_sc(X, Y,q,name='K',xlen = 7, yadjust = 0, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    peak = np.array([1,2,3,4,5,6,8,9,10,11,12,13,14,16])
    peaks = np.array([1,2,3,4,5,6,8,9,10,11,12,13,14,16])
    x = peaks**0.5*q
    y = np.array([ 0.8,0.8,0.8,0.7,0.6,0.5,0.4,0.3,0.17,0.17,0.17,0.17,0.1])
    y = np.array([ 0.8,0.8,0.8,0.7,0.6,0.5,0.4,0.3,0.7,0.6,0.5,0.4,0.3,0.15,0.1])
    for i in range(len(x)):
        #print x[i]
        #print y[i]
        plt.axvline(x=x[i],ymin=0,ymax=y[i],ls='--',color='grey')
        #print peak[i]
        #plt.text(x[i],y[i]+0.15,'$\sqrt{%i}$'%(peak[i]),size = 7,ha='center')

    x = np.array(X)
    y = np.array(Y)
    ax1.plot(x, y,color='#000080', linestyle='-', label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #Debye-Waller
    #plt.xlim(0,1.6)
    #plt.ylim(0,1.2)
    #annotate
    save_fig(save=save,fig_size=fig_size)
def plot_sf(X, Y,sf, name='K',xlen = 7, yadjust = 0, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    xx = np.arange(0,2,0.05)
    yy = np.zeros(xx.shape)
    from math import exp
    r=2.3
    #for i in range(xx.shape[0]):
    #    yy[i] = exp(-(r**2*xx[i]**2)/3)
    #ax1.plot(xx,yy,'--',color='green')
    x= [sf*i**0.5 for i in range(1,xlen+1)]
    #y = np.array([ 0.8,0.8,0.8,0.7,0.6,0.5,0.4,0.3,0.17])+0.2
    y = [0.85 for i in range(xlen)]
    #for i in range(xlen):
    #    print x[i]
    #    print y[i]
    #    plt.axvline(x=x[i],ymin=0,ymax=y[i],ls='--',color='grey')
    #    if i == 0:
    #        plt.text(x[i],y[i]+0.25,'q*',size = 7,ha='center')
    #    if i == 6 or i == 14:
    #        plt.text(x[i],y[i]+0.25,'$\sqrt{%i}$'%(i+1),size = 7,ha='center')

    #plt.text(1.45,1,name+'.',size=12,ha='center')
    x = np.array(X)
    y = np.array(Y)
    ax1.plot(x, y,color='#000080', marker='x', linestyle='', label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim(0,2)
    plt.ylim(0,1)
    #Debye-Waller
    #annotate
    save_fig(save=save,fig_size=fig_size)
