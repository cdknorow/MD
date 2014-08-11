import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import powerexpansioneos as peos
reload(peos)

#set the fonts
matplotlib.rc('font', **{'size':'14'})
counter = 0
colors = ['r','g','b','k','y','c','m']
choose = range(0,20,3)

def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def radial_plots(ax1):
    order = []
    os.chdir('histogram')
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file.split('.')[0][4:]),hist_file])
    order = sorted(order)
    count = 0
    order.append([order[-1][0]+1,'../hist_target.inp'])
    print 'numer of files =',len(order)
    for hist_file in order:
        if count in choose or hist_file == order[-1]:
            fid = open(hist_file[1],'r')
            dist = []
            r = []
            #Get Distances
            for line in fid.readlines():
                if line[0] != '#':
                    r.append(float(line.split()[0]))
                    dist.append(float(line.split()[1]))
            fid.close()
            if hist_file == order[-1]:
                ax1.plot(r,dist,lw='3',label=hist_file[1])
            else:
                ax1.plot(r,dist,label=hist_file[1])
           # ax1.set_xlim(1.4,3.5)
            #ax1.set_xlim(.8,2.2)
            ax1.set_ylim(0,4)
        count+=1
    os.chdir('../')

def full_radial_plots(ax1):
    order = []
    os.chdir('histogram_full')
    for hist_file in sorted(os.listdir('.')):
        order.append([float(hist_file.split('.dat')[0][9:]),hist_file])
    order = sorted(order)
    count = 0
    order.append([order[-1][0]+1,'../full_hist_target.inp'])
    print 'numer of files =',len(order)
    for hist_file in order:
        if count in choose or hist_file == order[-1]:
            fid = open(hist_file[1],'r')
            dist = []
            r = []
            #Get Distances
            for line in fid.readlines():
                dist.append(float(line.split()[1]))
                r.append(float(line.split()[0]))
            fid.close()
            if hist_file == order[-1]:
                ax1.plot(r,dist,lw='3',label=hist_file[1])
            else:
                ax1.plot(r,dist,label=hist_file[1])
        count+=1
    os.chdir('../')

def potential_plots(ax1):
    os.chdir('potential')
    order = []
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file.split('.')[0][9:]),hist_file])
    order = sorted(order)
    count = 0
    for hist_file in order:
        if count in choose or hist_file == order[-1]:
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
            ax1.plot(r,V,label=hist_file[1])
            #ax1.plot(r,F,label=hist_file)
            #ax1.set_xlim(.8,2.5)
            ax1.set_xlim(1.4,3.2)
            ax1.set_ylim(-.2,6)
        count+=1
    os.chdir('../')

if __name__ == "__main__":
    ax1 = get_plot()
    radial_plots(ax1)
    plt.legend(loc='upper right')
    plt.show()
    if True:
        #ax1 = get_plot()
        #full_radial_plots(ax1)
        #plt.legend(loc='upper right')
        #plt.show()
        pass
    else:
        ax1 = get_plot()
        full_radial_plots(ax1)
        plt.legend(loc='upper left')
        plt.show()
        ax1 = get_plot()
        potential_plots(ax1)
        plt.show()


