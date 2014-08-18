import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#set the fonts
matplotlib.rc('font', **{'size':'20'})
plt.close()
plt.figure(figsize=(20,12), dpi=100)
counter = 0
colors = ['r','g','b','k','y','c','m']

def pressure(ax1, title='pressure', xmin=0, xmax=1, index=0):
    global counter
    global colors
    file1 = 'pressure.log'
    fid = open(file1,'r')
    pressure = [[] for i in range(7)]
    fid.readline()
    Ptotal = []
    count = 0
    avg_count = 0
    delta = 5.0
    pavg = [0.0]
    for P in fid.readlines():
        if (count < xmax) and (count > xmin):
            for i in range(7):
                pressure[i].append(float(P.split()[i+1]))
            Ptotal.append((pressure[1][-1]+pressure[2][-1]+pressure[3][-1])/3.)
            pavg[-1] = pavg[-1] + pressure[0][-1]/delta
            avg_count+=1
            if avg_count > 5:
                avg_count = 0
                pavg.append(0.0)

        count+=1
    fid.close()

    del pavg[-1]
    #plot the pressures vs time
    if index == 0 :
        x = np.arange(len(pavg))*5
        ax1.plot(x, pavg, '%c'%colors[counter], label=title)
        counter+=1
    if index == 1 :
        x = range(len(pressure[1]))
        ax1.plot(x, pressure[0], '%c'%colors[counter], label=title)
        counter+=1
    if index == 2:
        x = range(len(pressure[1]))
        ax1.plot(x, pressure[0], label=title)
        ax1.plot(x, pressure[4], label='xy')
        ax1.plot(x, pressure[5], label='yz')
        ax1.plot(x, pressure[6], label='xz')
    #plot the pressures vs time
    if index == 3:
        x = range(len(pressure[1]))
        ax1.plot(x[5:], pressure[0][5:], label=title)
        ax1.plot(x[5:], pressure[1][5:], label='xx')
        ax1.plot(x[5:], pressure[2][5:], label='yy')
        ax1.plot(x[5:], pressure[3][5:], label='zz')
        ax1.plot(x[5:], Ptotal[5:], label='psum/3')


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
    if len(arg) > 1:
        count = 0
        if len(arg) == 2:
            fig1 = plt.figure(figsize=(16,8), dpi=100)
            x=1
            y=2
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            ax1 = fig1.add_subplot('%i%i%i'%(x,y,count))
            count+=1
            title = f.split('_')[0] +' '+ f.split('_')[-1]
            pressure(ax1,title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    index = int(config['-i']))
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
            pressure(ax1,title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    index = int(config['-i']))
            os.chdir('../')
    plt.legend(loc=2)
    plt.show()

