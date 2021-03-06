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
choose = [nruns/hist_num-1]
#hoose = range(20)

def Smooth(M):
    #Find the max
    for i in range(5):
        M_new = []
        for i in range(len(M)-1):
            M_new.append((M[i-1]+M[i]+M[i+1])/3.)
        M_new.append(M[-1])
        M = M_new
    return np.array(M_new)

def get_plot():
    fig1 = plt.figure(figsize=(20,8), dpi=100)
    ax1 = fig1.add_subplot(121)
    ax2 = fig1.add_subplot(122)
    return ax1,ax2

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
                    else:
                        if line[1] != 'r':
                            print 'L=',line.split()[3]
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
                ax1.set_ylim(0,np.max(dist)+.5)

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
                if big:
                    fid.readline()
                    fid.readline()
                for line in fid.readlines():
                    dist.append(float(line.split()[1]))
                    r.append(float(line.split()[0]))
                fid.close()
                dist = dist[10:]
                r = r[10:]
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
    for hist_file in sorted(os.listdir('../')):
        if hist_file[:16] == 'full_hist_target':
            order.append([counter,'../'+hist_file])
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
            ax1.plot(r,Smooth(V),'x-',label=hist_file[1])
            #ax1.plot(r,-np.gradient(V,r[1]-r[0]),'x-',label=hist_file[1]+'grad')
            ax1.set_xlim(r[0],r[-1])
        count+=1
    os.chdir('../')
    #ndna 200 sp3
    # epsilon = 1.0 
    # sigma = 10.5+3.13
    # p=14.5
    # l_c = 1.0
    #ndna 50 sp3
    # epsilon = 1.0
    # sigma = 10.5+1.82
    # p=10.5
    # l_c = 1.15
    #ndna 200 sp8
    #epsilon = 1.0
    #sigma = 10.5+8.75
    #p=11.75
    #l_c = 1.04
    # epsilon = 1.0
    # sigma = 10.5+8.5
    # p=11.75
    # l_c = 1.00
    #Grafted Polymer sp=4
    # epsilon = 1.0 
    # sigma = 5+3.7
    # p=8.5
    # l_c = 1.275
    #Grafted Polymer sp=8
    # epsilon = 1.0 
    # sigma = 5+7.
    # p=10.75
    # l_c = 1.42
    #Grafted Polymer sp=12
    epsilon = 160.0 
    sigma = 17.
    p=10
    n=10
    l_c = 1.55

    V=[]
    for r_x in r:
         #V.append(  epsilon * ( 1 / (np.sinh((r_x/sigma))**p)))
         #V.append(epsilon*(sigma/(r_x-sigma/l_c))*np.exp(-(r_x/sigma)**p))
         if r_x < sigma:
             V.append(epsilon*(1-r_x/sigma)**(5/2.))
         else:
             V.append(0)
         #V.append(epsilon*(sigma/(r_x))**n*np.exp(-(r_x/sigma)**p))
         #V.append(epsilon*(sigma/(r_x))**n + math.log(r_x/sigma))
         #V.append(-epsilon*5/18.*(30**(3/2.))*(math.log(r_x/sigma)))
         #V.append(epsilon*np.exp(-(r_x/sigma)**p))
         #V.append(epsilon*np.exp(-((r_x/sigma)**2)))
         #V.append(  epsilon * ( (sigma / r)**12))
    dt = 60
    #peos.rexp_Potential_fit(ax1, np.array(r[dt:-dt]), np.array(V[dt:-dt]),
    #        sigma=17,p=10,epsilon=1,l_c=1,show=True,umin=6)
    ax1.plot(r,V)
    print "Sigma:",sigma," p: ", p, " l_c: ", l_c
    plt.legend(loc=1)

if __name__ == "__main__":
    ax1,ax2 = get_plot()
    rmin=6
    rmax=20
    radial_plots(ax1)
    #plt.legend(loc='upper right')
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
        potential_plots(ax1)
        ax1.set_ylim(.01,5)
        potential_plots(ax2)
        ax2.set_yscale('log')
        #ax1.set_ylim(.01,200)
        ax2.set_ylim(.01,50)
    plt.show()


