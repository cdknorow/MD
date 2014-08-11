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
crystal = ['bcc','hcp','fcc','liq','mrco', 'solid']
############################################
def crystals(ax1, xmax=-1, title='',index = 3,xmin =-1):
    global counter
    global colors
    try:
        fid = open('crystal_count.dat','r')
        delta = int(fid.readline().split()[3])
        # bcc, hcp, fcc, liq, mirco
        M = [[],[],[],[],[],[]]
        for line in fid.readlines():
            total = float(int(line.split()[1]) + int(line.split()[2]) +
                    int(line.split()[3]) + int(line.split()[4]) +
                    int(line.split()[5]))
            M[0].append(int(line.split()[1])/total)
            M[1].append(int(line.split()[2])/total)
            M[2].append(int(line.split()[3])/total)
            M[3].append(int(line.split()[4])/total)
            M[4].append(int(line.split()[5])/total)
            M[5].append((M[0][-1]+M[1][-1]+M[2][-1]))
        x = range(0,len(M[0])*delta,delta)
    except:
        print 'no pickles'
        return 0
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
    ax1 = fig1.add_subplot(111)
    if len(arg) > 1:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] +' '+ f.split('_')[-1]
            crystals(ax1,title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    index = int(config['-i']))
            os.chdir('../')
    else:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] +' '+ f.split('_')[-1]
            for index in range(5):
                crystals(ax1, title = crystal[index], xmin = int(config['-x']), xmax = int(config['-r']),
                        index = index)
            os.chdir('../')
    plt.legend(loc=2)
    plt.show()
