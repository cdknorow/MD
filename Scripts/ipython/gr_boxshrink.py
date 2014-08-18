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
choose = range(1,20,2)

def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def radial_plots(ax1):
    os.chdir('histogram')
    #graph the plots in order
    def make_plots(order,choose=choose,big=False):
        count = 0
        for hist_file in order:
            if count in choose:
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
                    ax1.plot(r,dist,'-'+colors[count],lw='2',label=hist_file[1][14:-3])
                else:
                    ax1.plot(r,dist,'x'+colors[count],label=hist_file[1])
                ax1.set_ylim(0,4)
            count+=1
    order = []
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file.split('_')[0][4:]),hist_file])
    make_plots(sorted(order))
    print 'numer of files =',len(order)
    order = []
    counter = 0
    for hist_file in sorted(os.listdir('../')):
        if hist_file[:11] == 'hist_target':
            order.append([int(hist_file.split('.')[0][11:]),'../'+hist_file])
            counter+=1
    make_plots(sorted(order),big=True)
    ax1.set_xlim(rmin,rmax)
    ax1.set_ylim(0,3.5)
    os.chdir('../')

def full_radial_plots(ax1):
    os.chdir('histogram_full')
    def make_plots(order,choose=choose,big=False):
        count = 0
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
                if big:
                    ax1.plot(r,dist,colors[count],lw='1.5',label=hist_file[1][16:])
                else:
                    ax1.plot(r,dist,'x'+colors[count],label=hist_file[1])
            count+=1
    order = []
    for hist_file in sorted(os.listdir('.')):
        order.append([int(hist_file.split('_')[1][4:]),hist_file])
    order = sorted(order)
    print 'numer of files =',len(order)
    make_plots(order)
    order = []
    counter = 0
    for hist_file in sorted(os.listdir('../../')):
        if hist_file[:16] == 'full_hist_target':
            order.append([int(hist_file.split('.')[0][16:]),'../../'+hist_file])
            counter+=1
    make_plots(sorted(order),big=True)
    os.chdir('../')
    ax1.set_xlim(rmin,rmax*3)
    ax1.set_ylim(0,3.5)

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
    rmin = 7
    rmax = 18
    radial_plots(ax1)
    #plt.legend(loc='upper right')
    plt.show()
    if True:
        ax1 = get_plot()
        full_radial_plots(ax1)
        plt.legend(loc='upper right')
        plt.show()
        pass
    else:
        ax1 = get_plot()
        full_radial_plots(ax1)
        plt.legend(loc='upper right')
        plt.show()
        ax1 = get_plot()
        potential_plots(ax1)
        plt.show()


