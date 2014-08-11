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
colors = ['y','g','r','b','k','c','m','y-','g-','r-']
counter = 0
############################################
def msd(xmax=-1, time_scale=1, title='',xfit =-1,D=True):
    global counter
    global colors
    try:
        M = pickle.load(open('msd.pkl','r'))
        x = (M[0][0]+1).tolist()
        msd = M[1][0].tolist()
        x.insert(0,0.0)
        msd.insert(0,0.0)
        x = np.array(x)
        msd = np.array(msd)
    except:
        print 'no pickles'
        return 0
    #msd_fit
    #fit = []
    #for i in x:
    #plt.xscale('log')
    #plt.yscale('log')
    #    fit.append((2*3*D*(i) - tau*(1-math.e**(-i/tau)))**0.5)
    plt.plot(x[:xmax],msd[:xmax],'%s'%colors[counter], lw=5, label = title)
    #plot the diffusion coeffiencents
    #if D == 'True':
    #    def diffusion_coef(x,D,c):
    #        return ((6*D*x))**0.5+c
    #    from scipy.optimize import curve_fit
    #    popt, pcov = curve_fit(diffusion_coef, x[:xfit], msd[:xfit])
    #    plt.plot(x[:xmax],diffusion_coef(x[:xmax],popt[0],popt[1]),
    #            colors[counter], label='D=%.2f'%popt[0])
    if D == 'True':
        def diffusion_coef(x,D):
            return ((6*D*x))**0.5
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(diffusion_coef, x[:xfit], msd[:xfit])
        plt.plot(x[:xmax],diffusion_coef(x[:xmax],popt[0]),
                colors[counter], label='D=%.2f'%popt[0])
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
        msd(title = title, xmax = int(config['-r']), xfit = int(config['-x']))
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
                msd(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'])
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
                msd(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'])
                os.chdir('../../')
    else:
        msd(xmax = int(config['-r']), xfit = int(config['-x']))
    plt.xlabel = 'Time'
    plt.ylabel = 'MSD'
    plt.legend(loc=4)
    plt.show()

