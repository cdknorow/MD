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

def get_plot():
    fig1 = plt.figure(figsize=(20,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

############################################
def function(ax1, xmax=-1, time_scale=1, title='', binsize=.05):
    global counter
    global colors
    f = filter(lambda x: 'gaussian' in x, os.listdir('.'))
    index = []
    for i in f:
        index.append(i.split('.')[0].split('n')[1])
    index.sort(key=float)
    f_list = []
    for i in index:
        f_list.extend(filter(lambda x: x.split('.')[0].split('n')[1] == i, f))
    print f_list
    for pik in f_list:
        print pik
        M = pickle.load(open(pik,'r'))
        def shift(M,x,move,not_r=True):
            binsize = x[1]-x[0]
            s = np.argmax(M)
            delta = s-len(x)/2.
            #shift x array by delta
            if not_r:
                return np.where(delta < 0, x+binsize*abs(delta), x-binsize*delta)+move
            else:
                return x+move


        M[0][0] = shift(M[0][1],M[0][0],counter*.2)
        M[1][0] = shift(M[1][1],M[1][0],counter*.2)
        M[2][0] = shift(M[2][1],M[2][0],counter*.2)
        M[3][0] = shift(M[3][1],M[3][0],counter*.2,not_r=False)


        #expectation value of x
        #def expectation(f,x,binsize):
        #    expx = 0
        #    expx2 = 0
        #    for i in range(f.shape[0]):
        #        expx += x[i]*f[i]**2*binsize
        #    for i in range(f.shape[0]):
        #        expx2 += x[i]**2*f[i]**2*binsize
        #    return expx, expx2, expx2 - expx**2

        #print 'expecation value <x>, <x^2>, sigma^2'
        #print 'x', expectation(M[0][0],M[0][1],binsize)
        #print 'y', expectation(M[1][0],M[1][1],binsize)
        #print 'z', expectation(M[2][0],M[2][1],binsize)

        ax1.plot(M[0][0][:xmax],M[0][1][:xmax],'r', lw=1.5, label = 'x')
        ax1.plot(M[1][0][:xmax],M[1][1][:xmax],'g', lw=1.5, label = 'y')
        ax1.plot(M[2][0][:xmax],M[2][1][:xmax],'b', lw=1.5, label = 'z')
        ax1.plot(M[3][0][:xmax],M[3][1][:xmax],'k', lw=1.5, label = 'r')
        #ax1.set_xlim((M[0][1][-xmax],M[0][1][xmax]))
        counter+=1
    ax1.set_xlim((-.1,.1+(counter*.15)))
    #ax1.set_ylim((0,.06))
    ax1.set_xlabel(r'$\sigma$')
    ax1.set_ylabel('count')

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    ax1 = get_plot()
    config = {'-x':'-1','-r':'-1','-i':0}
    options, arg =  getopt.getopt(sys.argv[1:],'x:r:D:i:')
    config.update( dict(options))
    #add_subplot(x,x,1) # x rows, x columns, first plot
    function(ax1,title = '', time_scale = int(config['-x']), xmax = int(config['-r']),
            binsize = float(config['-i']))

    plt.show()

