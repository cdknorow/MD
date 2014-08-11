import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import powerexpansioneos as peos
reload(peos)

#set the fonts
matplotlib.rc('font', **{'size':'10'})
counter = 0
colors = ['r','g','b','k','y','c','m']
colors.extend(colors)
colors.extend(colors)
hist_num = 0
for i in os.listdir('.'):
    if i[:5] == 'hist_':
        hist_num+=1
nruns = len(os.listdir('histogram/'))
choose = [0,nruns/hist_num-1]

def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def radial_plots(ax1):
    os.chdir('histogram')
    #graph the plots in order
    def make_plots(order,choose=choose,big=False):
        counter = 0
        for hist_file in order:
            count = int(hist_file[1][-5])
            if hist_file[0] in choose or big:
                fid = open(hist_file[1],'r')
                dist = []
                r = []
                #Get Distances
                for line in fid.readlines():
                    if line[0] != '#':
                        r.append(float(line.split()[0]))
                        dist.append(float(line.split()[1]))
                fid.close()
                if big:
                    ax1.plot(r,dist,colors[counter],lw='3',label=hist_file[1][14:])
                    counter+=1
                else:
                    if int(hist_file[1][4]) == 0:
                        ax1.plot(r,dist,'o'+colors[count],label=hist_file[1])
                    else:
                        ax1.plot(r,dist,'v'+colors[count],label=hist_file[1])
                ax1.set_xlim(r[0],r[-1])

    order = []
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file[4:].split('_')[0]),hist_file])
    make_plots(sorted(order))
    print 'numer of files =',len(order)
    order = []
    counter = 0
    for hist_file in sorted(os.listdir('../')):
        if hist_file[:11] == 'hist_target':
            order.append([counter,'../'+hist_file])
            counter+=1
    print order
    make_plots(sorted(order),big=True)
    os.chdir('../')

def full_radial_plots(ax1,big=False):
    order = []
    os.chdir('histogram_full')
    print 'numer of files =',len(order)
    def make_plots(order,choose=choose,big=False):
        for hist_file in order:
            count = int(hist_file[1][-5])
            if hist_file[0] in choose or big:
                fid = open(hist_file[1],'r')
                dist = []
                r = []
                #Get Distances
                for line in fid.readlines():
                    dist.append(float(line.split()[1]))
                    r.append(float(line.split()[0]))
                fid.close()
                if big:
                    ax1.plot(r,dist,colors[count],lw='3',label=hist_file[1])
                else:
                    ax1.plot(r,dist,colors[count],label=hist_file[1])
                ax1.set_xlim(r[0],r[-1])
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file[9:].split('_')[0]),hist_file])
    order = sorted(order)
    make_plots(order)
    order = []
    counter = 0
    for hist_file in sorted(os.listdir('../../')):
        if hist_file[:16] == 'full_hist_target':
            order.append([counter,'../../'+hist_file])
            counter+=1
    make_plots(sorted(order),big=True)
    os.chdir('../')

def potential_plots(ax1):
    os.chdir('potential')
    order = []
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file.split('.')[0][9:]),hist_file])
    order = sorted(order)
    count = 0
    for hist_file in order:
        if count in choose:
            print count
            fid = open(hist_file[1],'r')
            F = []
            V = []
            r = []
            #Get Distances
            fid.readline()
            for line in fid.readlines():
                V.append(float(line.split()[1]))
                F.append(float(line.split()[2]))
                r.append(float(line.split()[0]))
            fid.close()
            ax1.plot(r,V,'x-',label=hist_file[1])
            #ax1.plot(r,-np.gradient(V,r[1]-r[0]),'x-',label=hist_file[1]+'grad')
            ax1.set_xlim(r[0],15)
        count+=1
    os.chdir('../')
    #epsilon = 1.0 
    #sigma = 11.8
    #p=9.2
    #V=[]
    #for r in np.linspace(rmin,rmax,300):
    #    V.append(  epsilon * ( 1 / (np.sinh((r/sigma))**p)))
    plt.legend(loc=1)

if __name__ == "__main__":
    ax1 = get_plot()
    rmin=6
    rmax=20
    radial_plots(ax1)
    ax1.set_ylim(0,3.5)
    plt.legend(loc='upper right')
    plt.show()
    if False:
        ax1 = get_plot()
        full_radial_plots(ax1)
        ax1.set_xlim(rmin,rmax*3)
        ax1.set_ylim(0,3.5)
        plt.legend(loc='upper right')
        plt.show()
    else:
        #ax1 = get_plot()
        #full_radial_plots(ax1)
        #ax1.set_xlim(0,rmax*3)
        #ax1.set_ylim(0,3.5)
        #plt.legend(loc='upper right')
        #plt.show()
        ax1 = get_plot()
        potential_plots(ax1)
        ax1.set_yscale('log')
        #ax1.set_ylim(.01,200)
        ax1.set_ylim(.01,20)
        plt.show()


