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
def params_ipy(width=900,mult=1):
    from math import sqrt
    #Set Latex Parameters 
    fig_width_pt = width  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]
    mpl.rc('font', family='serif')
    params = {'backend': 'ps',
              'axes.labelsize': 32,
              'text.fontsize': 32,
              'legend.fontsize': 32,
              'xtick.labelsize': 32,
              'ytick.labelsize': 32,
              'text.usetex': False,
              'axes.linewidth':5,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
def params_ipy_square(width=900,mult=1):
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
              'axes.labelsize': 32,
              'text.fontsize': 32,
              'legend.fontsize': 32,
              'xtick.labelsize': 32,
              'ytick.labelsize': 32,
              'text.usetex': False,
              'axes.linewidth':5,
              'figure.figsize': fig_size}
    mpl.rcParams.update(params)
    return fig_size
#save current figure in certaint format
def save_fig(save='fig',savefmt='png',fig_size=[3.251,2.00]):
    plt.savefig(save+'.'+savefmt,figsize=fig_size,dpi=300,
            bbox_inches='tight',format=savefmt,orientation='portrait',
            pad_inches=0.03)
