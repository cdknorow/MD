import sys
import os
import pickle
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import ipyplot
reload(ipyplot)


#set the fonts
matplotlib.rc('font', **{'size':'32'})
plt.close()
fig_size = ipyplot.params_ipy_square()
plt.figure(figsize=fig_size, dpi=100)
counter = 0
colors = ['y','g','r','b','k','c','--y','--g','--r','--b','--k']
crystal = ['bcc','hcp','fcc','other']

############################################
def crystals(ax1,xmin=-1, xmax=-1, title='', index = 3):
    global counter
    global colors
    fid = open('cna_crystal_count.dat','r')
    delta = int(fid.readline().split()[3])
    # bcc, hcp, fcc, other
    M = [[],[],[],[]]
    for line in fid.readlines():
        total = float(int(line.split()[1]) + int(line.split()[2]) +
                int(line.split()[3]) + int(line.split()[4]))
        M[0].append(int(line.split()[1])/total)
        M[1].append(int(line.split()[2])/total)
        M[2].append(int(line.split()[3])/total)
        M[3].append(int(line.split()[4])/total)
    x = range(0,len(M[0])*delta,delta)
    ax1.plot(x[:xmax],M[index][:xmax],'%s'%colors[counter], lw='5', label = title)
    ax1.set_ylim(-.05,1.1)
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
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    if len(arg) > 1:
        if len(arg) == 2:
            fig1 = plt.figure(figsize=(16,6), dpi=100)
            x=1
            y=2
            count = 1
        if len(arg) == 3 or len(arg) == 4:
            fig1 = plt.figure(figsize=(20,7), dpi=100)
            x=1
            y=4
            count = 3
        if len(arg) == 5 or len(arg) == 6:
            fig1 = plt.figure(figsize=(12,8), dpi=100)
            x=3
            y=2
            count = 5
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            ax1 = fig1.add_subplot('%i%i%i'%(x,y,count))
            count-=1
            counter = 0
            title = f.split('_')[9] 
            ax1.set_title(r'$\rho = $'+title)
            for index in range(4):
                crystals(ax1, title = crystal[index], xmin = int(config['-x']), xmax = int(config['-r']),
                        index = index)
            os.chdir('../')
    else:
        ax1 = fig1.add_subplot(111)
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] +' '+ f.split('_')[-1]
            for index in range(4):
                crystals(ax1, title = crystal[index], xmin = int(config['-x']), xmax = int(config['-r']),
                        index = index, text=title)
            os.chdir('../')
    plt.legend(loc=6)
    plt.show()



