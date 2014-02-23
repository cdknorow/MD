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
plt.figure(figsize=(20,12), dpi=100)
counter = 0
colors = ['r','g','b','k','y','c','m']
############################################
def self_int(xmax=-1, time_scale=1, title='',xfit =-1,D=True):
    global counter
    global colors
    try:
        M = pickle.load(open('si_scat.pkl','r'))
        x = np.array(M[0])
        si = M[1]
    except:
        print 'no pickles'
        return 0
    plt.plot(x[:xmax],si[:xmax],'%cx'%colors[counter], label = title)
    plt.yscale('log')
    #plot the diffusion coeffiencents
    if D == 'True':
        from scipy.optimize import curve_fit
        def stretched_exp(x,B,ta):
            return .99 * np.exp(-(x/ta)**B)
        try:
            popt, pcov = curve_fit(stretched_exp, x[:xfit], si[:xfit])
            plt.plot(x[:xmax],stretched_exp(x[:xmax],popt[0],popt[1]),
                    colors[counter], label='B=%.2f'%popt[0])
        except:
            'fit failed'
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
        self_int(title = title, xmax = int(config['-r']), xfit = int(config['-x']))
        os.chdir('../')
    elif len(arg) > 0:
        if config['-l']=='single':
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    dirs.append(i)
            for i in arg:
                f = dirs[int(i)]
                print f
                os.chdir(f)
                if f[0] =='F':
                    title = 'rho = %s'%f.split('_')[9]
                if f[0] =='s':
                    title = 'rho = %s'%f.split('_')[10]
                #if f[0] =='F':
                #    title = 'N = %s'%f.split('_')[13]
                #if f[0] =='s':
                #    title = 'N = %s'%f.split('_')[14]
                self_int(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'])
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
                print r
                if r[0] =='F':
                    title = 'rho = %s'%r.split('_')[9]
                if r[0] =='s':
                    title = 'rho = %s'%r.split('_')[10]
                if r[0] =='F':
                    title = 'N = %s'%r.split('_')[13]
                if r[0] =='s':
                    title = 'N = %s'%r.split('_')[14]
                self_int(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'])
                os.chdir('../../')
    else:
        self_int(xmax = self_int(config['-r']), xfit = int(config['-x']))
    plt.xlabel = 'Time'
    plt.ylabel = 'MSD'
    plt.legend(loc=3)
    plt.show()


