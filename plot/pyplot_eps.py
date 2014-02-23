## \package MD.analysis.pyplot 
# \brief This module is used for making simple plots on the fly using matplotlib

import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as plines

def color_pick(name):
    if name == 'green':
        return '#325431'
    if name == 'orange':
        return '#db6214'
    if name == 'blue':
        return '#265787'
    if name == 'red':
        return '#850000'
    if name == 'burnt':
        return '#5e0d0d'
    if name == 'brown':
        return '#5d2e09'
    if name == 'purple':
        return '#543d61'
def interp(bins,hist):
    from scipy import interpolate
    s = interpolate.InterpolatedUnivariateSpline(bins, hist)
    xs = np.arange(bins[0],bins[-1],bins[0]-bins[1]/10)
    hist_s=s(xs)
    return xs, hist_s
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
              'text.fontsize': 9,
              'legend.fontsize': 6,
              'xtick.labelsize': 7,
              'ytick.labelsize': 7,
              'text.usetex': False,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
#save current figure in certaint format
def save_fig(save='fig',savefmt='eps',fig_size=[3.251,2.00]):
    plt.savefig(save+'.'+savefmt,figsize=fig_size,dpi=300,
            bbox_inches='tight',format=savefmt,orientation='portrait',
            pad_inches=0.03)
#### The Main Plotting function
def plot(X, Y, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, y,color=color_pick('blue'), marker=marker, linestyle=linestyle, label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim()
    #plt.xlim((0,1e8))
    if showleg:
        ax1.legend(loc=2)
    if lim != False:
        plt.xlim(lim)
        plt.ylim([-1,100])
    save_fig(save=save,fig_size=fig_size)
def plot2(X1,Y1,X2,Y2,xlabel='x',ylabel='y',label1='',
        label2='',save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, 's', color = '#990000',linewidth = 1, label=label1)
    ax1.plot(X2, Y2,'--',color = 'k',linewidth = 1, label=label2)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if showleg:
        ax1.legend(loc=2)
    if lim != False:
        plt.xlim(lim)
        plt.ylim(0,1)

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
    #plt.ylim([0,1.2])
    if lim != False:
        plt.xlim(lim)
    if showleg:
        ax1.legend(loc=1)
    save_fig(save=save,fig_size=fig_size)
def plot4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4='',save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000', label=label1)
    ax1.plot(X2, Y2,color = '#CC9900', label=label2)
    ax1.plot(X3, Y3, color = '#0099FF', label=label3)
    ax1.plot(X4, Y4, color = '#006666', label=label4)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim([0,1.8e8])
    plt.ylim([0,1])
    if showleg:
        ax1.legend(loc=4)
    save_fig(save=save,fig_size=fig_size)
def plotextra(X1,Y1,X2,Y2,X3,Y3,X4,Y4,x1extra,y1extra,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4='',save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000', label=label1)
    ax1.plot(x1extra, y1extra,'x', color = '#990000')
    ax1.plot(X2, Y2 ,color = '#CC9900', label=label2)
    ax1.plot(X3, Y3, color = '#0099FF', label=label3)
    ax1.plot(X4, Y4, color = '#006666', label=label4)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim([0,2.2e8])
    plt.ylim([0,1])
    if showleg:
        ax1.legend(loc=4)
    save_fig(save=save,fig_size=fig_size)
def plot_multi(X, Y, label, D,step, xlabel='x', ylabel='y',save='fig'):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for i in range(y.shape[0]):
        ax1.plot(x, y[i], label = label[i])
    for i,d in enumerate(D):
        fit = (6*d*X)**0.5+step[i]
        ax1.plot(X, fit)
    ax1.legend(loc=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    save_fig(save=save,fig_size=fig_size)
#### histogram plottiong
def plot_bar(X,Y,x2,y2,xlabel='x',ylabel='y',title='',label='',
        save='fig1', showleg=False, savefmt='png', width=0.4):
    #Plot Data
    fig_size=params()
    y_array = np.array(Y)
    x_array = np.array(X)
    fig1=plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.bar(x_array,y_array,width)
    ax1.plot(x2,y2,color='k')
    plt.axvline(x=4.4,ymin=0,ymax=35,ls='--',color='r')
    ax1.text(7,70,r'$\lambda^2=4$',size=8)
    plt.xlim([1,15])
    if showleg:
        ax1.legend(loc=2)
    ax1.set_xlabel(xlabel,labelpad=2)
    ax1.set_ylabel(ylabel,labelpad=2)
    ax1.axes.xaxis.get_major_locator()._nbins=7
    ax1.axes.yaxis.get_major_locator()._nbins=5
    save_fig(fig_size=fig_size,save=save,savefmt='eps')
##################### paper specific
def plot_crystal_time(X1,Y1,X2,Y2,X3,Y3,X4,Y4,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4='',save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000', label=label1)
    ax1.plot(X2, Y2,color = '#CC9900', label=label2)
    ax1.plot(X3, Y3, color = '#0099FF', label=label3)
    ax1.plot(X4, Y4, color = '#006666', label=label4)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim([0,1.8e8])
    plt.ylim([0,1])
    if showleg:
        ax1.legend(loc=4)
    save_fig(save=save,fig_size=fig_size)
def plot_num_crystal(X1,Y1,X2,Y2,X3,Y3,xlabel='x',ylabel='y',label1='',
        label2='',label3='',save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(311)
    ax2 = fig1.add_subplot(312)
    ax3 = fig1.add_subplot(313)
    ax1.plot(X1, Y1, color = color_pick('blue'), label=label1)
    plt.ylim([0,450])
    plt.xlim([0,2e8])
    plt.ylim([0,450])
    plt.xlim([0,2e8])
    #ax2.fill_between(X2,0, Y2, color = '#660000', facecolor = 'red', linewidth =0.25, label=label2)
    ax2.plot(X2, Y2, color = color_pick('blue'), label=label2)
    plt.ylim([0,25])
    ax3.plot(X3,Y3, color_pick('red'),label=label3)
    plt.ylim([0,10])
    #set the labels to be shown
    ax3.set_xlabel(xlabel)
    ax2.set_ylabel('NPs')
    ax1.set_ylabel('NPs')
    ax3.set_ylabel('Domains')
    #set the axis visibl
    ax2.xaxis.set_visible(False)
    ax1.xaxis.set_visible(False)
    #adjust the spacing
    plt.subplots_adjust(wspace=0.1)
    plt.subplots_adjust(hspace=0.2)
    #adjust the number of ticks
    ax1.axes.yaxis.get_major_locator()._nbins=5
    ax3.axes.yaxis.get_major_locator()._nbins=4
    ax2.axes.yaxis.get_major_locator()._nbins=4
    #set the legend location
    #ax1.legend(loc=2)
    #ax2.legend(loc=1)
    #ax3.legend(loc=1)
    ax1.annotate('Critical Nucleus',xy=(1e8,75),xytext=(1e8,250), size = 7, 
            ha = 'right', arrowprops = dict(facecolor =
        'black',shrink=0.1, width = 0.25, headwidth = 2.0))
    ax1.text(1.31e8,200,'Largest Crystal',size=8)
    ax2.text(1.3e8,15,'Second Largest',size=8)
    ax3.text(1.15e8,6,'Number of Domains',size=8)

    save_fig(save=save,fig_size=fig_size)
def plot_num_crystal4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,xlabel='x',ylabel='y',label1='',
        label2='',label3='',label4='',save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params(mult=1.4)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(411)
    ax2 = fig1.add_subplot(412)
    ax3 = fig1.add_subplot(413)
    ax4 = fig1.add_subplot(414)
    ax1.plot(X1, Y1, color = color_pick('blue'), label=label1)
    plt.ylim([0,450])
    plt.xlim([0,2e8])
    #ax2.fill_between(X2,0, Y2, color = '#660000', facecolor = 'red', linewidth =0.25, label=label2)
    ax2.plot(X2, Y2, color = color_pick('blue'), label=label2)
    plt.ylim([0,25])
    ax3.plot(X3,Y3, color_pick('red'),label=label3)
    plt.ylim([0,10])
    ax4.plot(X4, Y4, color = color_pick('orange'), lw=2, label=label4)
    Y4l=[]
    for i in range(X4.shape[0]):
        Y4l.append(8.57)
    ax4.plot(X4, Y4, color = color_pick('orange'), lw=2, label=label4)
    ax4.plot(X4, Y4l,'--', color = 'k', lw=1, label=label4)
    plt.ylim([7.5,9])
    #set the labels to be shown
    ax3.set_xlabel(xlabel)
    ax2.set_ylabel('NPs')
    ax1.set_ylabel('NPs')
    ax3.set_ylabel('Domains')
    ax4.set_ylabel('f$(h$)')
    #set the axis visibl
    ax2.xaxis.set_visible(False)
    ax1.xaxis.set_visible(False)
    ax3.xaxis.set_visible(False)
    #adjust the spacing
    plt.subplots_adjust(wspace=0.1)
    plt.subplots_adjust(hspace=0.2)
    #adjust the number of ticks
    ax1.axes.yaxis.get_major_locator()._nbins=4
    ax3.axes.yaxis.get_major_locator()._nbins=4
    ax2.axes.yaxis.get_major_locator()._nbins=4
    ax4.axes.yaxis.get_major_locator()._nbins=3
    #set the legend location
    #ax1.legend(loc=2)
    #ax2.legend(loc=1)
    #ax3.legend(loc=1)
    ax1.annotate('Critical Nucleus',xy=(1e8,75),xytext=(1e8,250), size = 7, 
            ha = 'right', arrowprops = dict(facecolor =
        'black',shrink=0.1, width = 0.25, headwidth = 2.0))
    ax1.text(1.2e8,200,'Largest Crystal',size=8)
    ax2.text(1.2e8,15,'Second Largest',size=8)
    ax3.text(1.05e8,6,'Number of Domains',size=8)
    ax4.text(0.9e8,4,'Hybridizations per NP',size=8)

    save_fig(save=save,fig_size=fig_size)
def plot_defects(X1,Y1,X2,Y2,X3,Y3,xlabel='x',ylabel='y',label1='',
        label2="",label3="",save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plt.xlim([0.5e8,2.1e8])
    ax1.plot(X2, Y2, color = color_pick('green'), label=label2)
    #ax1.plot(X2[0:-1:100], Y2[0:-1:100], 'x', ms=2.0,color = color_pick('green'))
    ax1.plot(X1, Y1, color = color_pick('orange'),lw=0.5, label = label1)
    #ax1.plot(X1[0:-1:200], Y1[0:-1:200],'|', color = color_pick('orange'),ms=4.0 )
    ax1.plot(X3, Y3, color = color_pick('blue'), label=label3)
    #plt.xlim([1e8,1.6e8])
    #plt.ylim([0,40])
    plt.xlim([0.5e8,1.3e8])
    plt.ylim([0,25])
    plt.axvline(x=1.3e8,ymin=0,ymax=25,ls='--', color='black')
    ax1.annotate('K',xy=(1.4e8,20),xytext=(1.5e8,22),size = 12, arrowprops = dict(facecolor =
        'black', shrink = 0.1, width = 0.2, headwidth = 2))
    ax1.annotate('L',xy=(1.5e8,13),xytext=(1.5e8,16),size = 12,arrowprops = dict(facecolor =
        'black', shrink = 0.1, width = 0.2, headwidth = 2),
        horizontalalignment = 'center')
    ax1.annotate('M',xy=(1.9e8,7),xytext=(1.9e8,10),size = 12, arrowprops = dict(facecolor =
        'black', shrink = 0.1, width = 0.2, headwidth = 2),horizontalalignment = 'center')
    ax1.annotate('',xy=(1.15e8,3.25),xytext=(1.23e8,3.25),
            size = 6, arrowprops = dict(facecolor =
        'black', shrink = 0.1, width = 0.25, headwidth = 1.25))
    plt.text(1.23e8,2.25,'Vacancy-Interstitial Pair',size = 6)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.legend(loc=6)
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
def plot_sf_sc(X, Y,sf, name='K',xlen = 7, yadjust = 0, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    xx = np.arange(0,2,0.05)
    yy = np.zeros(xx.shape)
    from math import exp
    r=2.3
    for i in range(xx.shape[0]):
        yy[i] = exp(-(r**2*xx[i]**2)/3)
    ax1.plot(xx,yy,'--',color='green')
    peak = np.array([1,2,3,4,5,6,8,9,10,11,12,13,14,16])
    peaks = np.array([1,2,3,4,5,6,8,9,10,11,12,13,14,16])
    x = peaks**0.5*sf
    y = np.array([ 0.8,0.8,0.8,0.7,0.6,0.5,0.4,0.3,0.17,0.17,0.17,0.17,0.1])
    y = np.array([ 0.8,0.8,0.8,0.7,0.6,0.5,0.4,0.3,0.7,0.6,0.5,0.4,0.3,0.15,0.1])
    for i in range(len(x)):
        print x[i]
        print y[i]
        plt.axvline(x=x[i],ymin=0,ymax=y[i],ls='--',color='grey')
        print peak[i]
        plt.text(x[i],y[i]+0.15,'$\sqrt{%i}$'%(peak[i]),size = 7,ha='center')

    plt.text(1.45,1,name+'.',size=12,ha='center')
    x = np.array(X)
    y = np.array(Y)
    ax1.plot(x, y,color='#000080', marker='x', linestyle='', label=label)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #Debye-Waller
    plt.xlim(0,1.6)
    plt.ylim(0,1.2)
    #annotate
    save_fig(save=save,fig_size=fig_size)
def plot_msd_int(X, Y, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, y,color=color_pick('burnt'), marker=marker, linestyle=linestyle, label=label)
    plt.axvline(x=1.3e8,ymin=0,ymax=100,ls='--', color='black')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim([0.5e8,2.1e8])
    save_fig(save=save,fig_size=fig_size)
def plot_sc(X, Y, linestyle='-', marker='', xlabel='x', ylabel='y',
        label='', save='fig', showleg=False, lim=False):
    #just to make sure they are arrays 
    x = np.array(X)
    y = np.array(Y)
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, y,color=color_pick('blue'), marker=marker, linestyle=linestyle, label=label)
    plt.axvline(x=1.3e8,ymin=0,ymax=1,ls='--', color='black')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if showleg:
        ax1.legend(loc=2)
    if lim != False:
        plt.xlim(lim)
        plt.ylim(0,1)
    save_fig(save=save,fig_size=fig_size)
def plot_defects_cycle(X1,Y1,X2,Y2,X3,Y3,X4,Y4,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4="",save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params(mult=1)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X2, Y2, color = color_pick('green'), label=label2)
    #ax1.plot(X2[0:-1:100], Y2[0:-1:100], 'x', ms=2.0,color = color_pick('green'))
    ax1.plot(X1, Y1, color = color_pick('orange'),lw=0.5, label = label1)
    ax1.plot(X1[0:-1:25], Y1[0:-1:25],'|', color = color_pick('orange'),
            ms = 4.0, mew = 1.5)
    ax1.plot(X3, Y3, color = color_pick('blue'), label=label3)
    ax1.plot(X4, Y4,'--', color = color_pick('purple'), label=label3)
    #plt.xlim([0.5e8,2.1e8])
    #plt.xlim([0.5e8,2.1e8])
    #plt.ylim([0,25])
    #plt.axvline(x=1.3e8,ymin=0,ymax=25,ls='--', color='black')
    #ax1.annotate('K',xy=(1.4e8,20),xytext=(1.5e8,22),size = 'small', arrowprops = dict(facecolor =
    #    'black', shrink = 0.2, width = 0.5, headwidth = 0.75))
    #ax1.annotate('L',xy=(1.5e8,13),xytext=(1.6e8,16),size = 'small', arrowprops = dict(facecolor =
    #    'black', shrink = 0.2, width = 0.5, headwidth = 0.75),
    #    horizontalalignment = 'center')
    #ax1.annotate('M',xy=(1.9e8,7),xytext=(1.9e8,10),size = 'small', arrowprops = dict(facecolor =
    #    'black', shrink = 0.2, width = 0.5, headwidth = 0.75),horizontalalignment = 'center')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.legend(loc=1)
    save_fig(save=save,fig_size=fig_size)
def plot_cycle_hybrid(X1,Y1,xlabel='x',ylabel='y',label1='',save='fig',showleg=False):
    #just to make sure they are arrays 
    fig_size=params(mult=0.4)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = color_pick('blue'), label = label1)
    #ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.xaxis.set_visible(False)
    #ax1.legend(loc=1)
    save_fig(save=save,fig_size=fig_size)
def plot_defects_cycle_all(X1,Y1,X2,Y2,X3,Y3,X4,Y4,X5,Y5,xlabel='x',ylabel='y',label1='',
        label2="",label3="",label4="",save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params(mult=1)
    fig1 = plt.figure()
    ax1 = plt.subplot2grid((4,3),(0,0), colspan=3)
    ax1.plot(X5, Y5, color = color_pick('burnt'), label=label3)
    ax1.xaxis.set_visible(False)
    ax1.set_ylabel('$T$')
    ax1.axes.yaxis.get_major_locator()._nbins=5
    ax2 = plt.subplot2grid((4,3),(1,0), rowspan=3, colspan =3)
    ax2.plot(X2, Y2, color = color_pick('green'), label=label2)
    ax2.plot(X1, Y1, color = color_pick('orange'),lw=0.5, label = label1)
    ax2.plot(X3, Y3, color = color_pick('blue'), label=label3)
    ax2.plot(X4, Y4,'--', color = 'k')
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    ax2.annotate('',xy=(2e7,6),xytext=(2e7,4),size = 'small', arrowprops = dict(facecolor =
        'black', shrink = 0.1, width = 0.25, headwidth = 2))
    ax2.text(1.8e7,3,'substitutionals \nat T=1.115',size=6)
    plt.subplots_adjust(wspace=0.2)
    ax2.legend(loc=6,frameon=False,labelspacing=0.25,borderpad=0.05)
    ax2.axes.yaxis.get_major_locator()._nbins=5
    save_fig(save=save,fig_size=fig_size)
def plot_neighbors(X1,Y1,X2,Y2,xlabel='x',ylabel='y',label1='',
        label2='',save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1, color = '#990000',linewidth = 1, label=label1)
    ax1.plot(X2, Y2, color = '#336622',linewidth = 1, label=label2)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.legend(loc=2)

    save_fig(save=save,fig_size=fig_size)
def plot_MSD(X1,Y1,X2,Y2,X3,Y3,D,xlabel='x',ylabel='y',label1='',
        label2="",label3="",save='fig',showleg=False,lim=False):
    #just to make sure they are arrays 
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(X1, Y1,'x', color = '#990000', label='a')
    ax1.plot(X2, Y2,'+', color = '#CC9900', label='b')
    ax1.plot(X3, Y3,'s', color = '#0099FF', label='c')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    step = [4,6,3]
    for i,d in enumerate(D):
        fit = (6*d*X1)**0.5+step[i]
        ax1.plot(X1, fit)
    ax1.legend(loc=2)
    save_fig(save=save,fig_size=fig_size)
