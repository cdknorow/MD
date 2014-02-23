import os
import sys
import math
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
#save current figure in certaint format
def save_fig(save='fig',savefmt='png',fig_size=[3.251,2.00]):
    plt.savefig(save+'.'+savefmt,figsize=fig_size,dpi=300,
            bbox_inches='tight',format=savefmt,orientation='portrait',
            pad_inches=0.03)
def plot(X,Y,barx,bary,ylabel='y',xlabel='',label = [],save='nnplot_bccfit',showleg=False,lim=False,
        width = 0.1):
    #just to make sure they are arrays 
    color = ['red','green','blue','orange']
    fig_size=params()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for i in range(len(X)):
        ax1.plot(X[i], Y[i], color = color[i], label=label[i])
    ax1.bar(barx,bary,width,color='black')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.xlim(0,40)
    plt.ylim([0,0.02])
    ax1.legend(loc=1)
    fig1.show()
    save_fig(save=save,fig_size=fig_size)
a = 11
q = 9.42
a = 2*q/3.**0.5
x = [a*3.**0.5/2,a,2.**0.5*a,3.**0.5/2*a,a*3**0.5]
print x
y = [1,1,1,1,1]
fid = open('nnbakos.txt','r')
M = fid.readlines()
R = [[] for i in range(5)]
for line in M:
    s = line.split()
    for i in range(len(R)):
        R[i].append(float(s[i]))

X = [R[0],R[0],R[0],R[0]]
Y = [R[1],R[2],R[3],R[4]]
plot(X,Y,x,y,label=['A-A','A-B','A-C','A-D'])

