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
def crystals(xmax=-1, time_scale=1, title='',xfit =-1,D=True,index = 3):
    global counter
    global colors
    crystal = ['bcc ','hcp ','fcc ','liq ','mirco ', 'solid']
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
    plt.plot(x[:xmax],M[index][:xmax],'%c'%colors[counter], label = crystal[index] + title)
    plt.ylim(-.05,1.1)
    counter+=1

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    config = {'-f':'','-x':'-1','-r':'-1','-D':'True','-l':'single','-i':3}
    options, arg =  getopt.getopt(sys.argv[1:],'f:x:r:D:l:i:')
    config.update( dict(options))
    if config['-f'] != '':
        os.chdir(config['-f'])
        title = config['-f'].split('_')[10]
        crystals(title = title, xmax = int(config['-r']), xfit = int(config['-x']))
        os.chdir('../')
    elif len(arg) > 1:
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
                crystals(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'], 
                        index = int(config['-i']))
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
                #if r[0] =='F':
                #    title = 'N = %s'%r.split('_')[13]
                #if r[0] =='s':
                #    title = 'N = %s'%r.split('_')[14]
                crystals(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'],
                        index = int(config['-i']))
                os.chdir('../../')
    else:
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
                for i in range(5):
                    crystals(title = title, xmax = int(config['-r']),
                            xfit = int(config['-x']), D = config['-D'], 
                            index = i)
                os.chdir('../')
    plt.xlabel = 'Time'
    plt.ylabel = 'MSD'
    plt.legend(loc=2)
    plt.show()


