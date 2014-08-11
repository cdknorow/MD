import sys
import os
import pickle
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#set the fonts
matplotlib.rc('font', **{'size':'20'})
plt.close()
plt.figure(figsize=(20,20), dpi=100)
colors = ['b','r','g','b','k','y','c','m']
counter = 0
############################################
def function(ax1, xmax=-1, time_scale=1, title='', binsize=.05):
    global counter
    global colors
    try:
        M = pickle.load(open('gaussian.pkl','r'))
    except:
        print 'no pickles'
        return 0
    def shift(M,x):
        s = M.argmax()
        delta = s-len(x)/2
        #shift x array by delta
        if delta < 0:
            x = x+binsize*abs(delta)
        if delta > 0:
            x = x-binsize*delta
        return x
    M[0][1] = shift(M[0][0],M[0][1])
    M[1][1] = shift(M[1][0],M[1][1])
    M[2][1] = shift(M[2][0],M[2][1])

    #expectation value of x
    def expectation(f,x,binsize):
        expx = 0
        expx2 = 0
        for i in range(f.shape[0]):
            expx += x[i]*f[i]**2*binsize
        for i in range(f.shape[0]):
            expx2 += x[i]**2*f[i]**2*binsize
        return expx, expx2, expx2 - expx**2

    #print 'expecation value <x>, <x^2>, sigma^2'
    #print 'x', expectation(M[0][0],M[0][1],binsize)
    #print 'y', expectation(M[1][0],M[1][1],binsize)
    #print 'z', expectation(M[2][0],M[2][1],binsize)

    ax1.plot(M[0][1][:xmax],M[0][0][:xmax],'%c'%colors[counter], lw=5, label = 'x')
    ax1.plot(M[1][1][:xmax],M[1][0][:xmax],'%c'%colors[counter+1], lw=5, label = 'y')
    ax1.plot(M[2][1][:xmax],M[2][0][:xmax],'%c'%colors[counter+2], lw=5, label = 'z')
    #ax1.set_xlim((M[0][1][-xmax],M[0][1][xmax]))
    ax1.set_xlim((-2,2))
    ax1.set_ylim((0,.06))
    ax1.set_xlabel(r'$\sigma$')
    ax1.set_ylabel('count')
    plt.legend(loc=2)
    counter+=1

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    plt.close()
    config = {'-x':'-1','-r':'-1','-i':0}
    options, arg =  getopt.getopt(sys.argv[1:],'x:r:D:i:')
    config.update( dict(options))
    #add_subplot(x,x,1) # x rows, x columns, first plot
    if len(arg) > 1:
        if len(arg) == 2:
            fig1 = plt.figure(figsize=(16,8), dpi=100)
            x=1
            y=2
        if len(arg) == 3 or len(arg) == 4:
            fig1 = plt.figure(figsize=(12,8), dpi=100)
            x=2
            y=2
        if len(arg) == 5 or len(arg) == 6:
            fig1 = plt.figure(figsize=(12,8), dpi=100)
            x=3
            y=2
        dirs = []
        count = 0
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[9]
            ax1 = fig1.add_subplot('%i%i%i'%(x,y,count))
            count+=1
            counter = 0
            ax1.set_title(r'$\rho$='+title)
            function(ax1,title = title, time_scale = int(config['-x']), xmax = int(config['-r']),
                    binsize = float(config['-i']))
            os.chdir('../')
    else:
        fig1 = plt.figure(figsize=(12,8), dpi=100)
        dirs = []
        ax1 = fig1.add_subplot(111)
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] +' '+ f.split('_')[-1]
            function(ax1,title = title, time_scale = int(config['-x']), xmax = int(config['-r']),
                    binsize = float(config['-i']))
            os.chdir('../')
    plt.show()

