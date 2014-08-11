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
def function(xmax=-1, time_scale=1, title='',xfit =-1,D=True, binsize=.05):
    global counter
    global colors
    try:
        M = pickle.load(open('radial.pkl','r'))
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
        print delta
        return x

    #expectation value of x
    def expectation(f,x,binsize):
        expx = 0
        expx2 = 0
        for i in range(f.shape[0]):
            expx += x[i]*f[i]**2*binsize
        for i in range(f.shape[0]):
            expx2 += x[i]**2*f[i]**2*binsize
        return expx, expx2, expx2 - expx**2


    plt.plot(M[0][1][:xmax],M[0][0][:xmax],'%c'%colors[counter], lw=5, label = 'x')
    counter+=1

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    config = {'-f':'','-x':'-1','-r':'-1','-D':'True','-l':'single'}
    options, arg =  getopt.getopt(sys.argv[1:],'f:x:r:D:l:')
    config.update( dict(options))
    if config['-f'] != '':
        os.chdir(config['-f'])
        title = config['-f'].split('_')[10]
        function(title = title, xmax = int(config['-r']), xfit = int(config['-x']))
        os.chdir('../')
    elif len(arg) > 0:
        if config['-l']=='single':
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    dirs.append(i)
            for i in arg:
                f = dirs[int(i)]
                os.chdir(f)
                if f[0] =='F':
                    title = 'rho = %s'%f.split('_')[9]
                if f[0] =='s':
                    title = 'rho = %s'%f.split('_')[10]
                #if f[0] =='F':
                #    title = 'N = %s'%f.split('_')[13]
                #if f[0] =='s':
                #    title = 'N = %s'%f.split('_')[14]
                function(title = title, binsize = float(config['-r']),
                        )
                print f
                os.chdir('../')
        else:
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    for j in sorted(os.listdir(i)):
                        if os.path.isdir(i+'/'+j):
                            dirs.append(i+'/'+j)
            for i in arg:
                f = dirs[int(i)]
                os.chdir(f)
                r = f.split('/')[1]
                if r[0] =='F':
                    title = r'$\rho$ = %s'%r.split('_')[9]
                if r[0] =='s':
                    title = r'$\rho$ = %s'%r.split('_')[10]
                if r[0] =='F':
                    title = 'N = %s'%r.split('_')[13]
                if r[0] =='s':
                    title = 'N = %s'%r.split('_')[14]
                print r
                function(title = title, xmax = int(config['-r']))
                os.chdir('../../')
    else:
        function(xmax = int(config['-r']))
    plt.xlabel = 'Sigma'
    plt.ylabel = 'Count'
    plt.legend(loc=4)
    plt.show()


