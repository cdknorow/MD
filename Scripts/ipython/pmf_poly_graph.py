import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import powerexpansioneos as peos
import h5py
reload(peos)

#set the fonts
matplotlib.rc('font', **{'size':'10'})
counter = 0
colors = ['r','g','b','k','y','c','m']
colors.extend(colors)
colors.extend(colors)

fig1 = plt.figure(figsize=(20,8), dpi=100)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

def Smooth(M):
    #Find the max
    for i in range(5):
        M_new = []
        for i in range(len(M)-1):
            M_new.append((M[i-1]+M[i]+M[i+1])/3.)
        M_new.append(M[-1])
        M = M_new
    return np.array(M_new)


def gen_plots():
    order = []
    f = filter(lambda x: '.pmf' in x, os.listdir('.'))
    count = 0
    for hist_file in sorted(f):
        fid = open(hist_file,'r')
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
        ax1.plot(r,Smooth(V),'x-',label=hist_file.split('.')[0])


    # epsilon = 1.0 
    # sigma = 16
    # p=14
    # n=12
    # l_c = 1.42
    # V=[]
    # for r_x in r:
    #      #V.append(  epsilon * ( 1 / (np.sinh((r_x/sigma))**p)))
    #      #V.append(epsilon*(sigma/(r_x-sigma/l_c))*np.exp(-(r_x/sigma)**p))
    #      #V.append(epsilon*(sigma/(r_x))**n*np.exp(-(r_x/sigma)**p))
    #      #V.append(epsilon*(sigma/(r_x))**n + math.log(r_x/sigma))
    #      V.append(-epsilon*5/18.*(30**(3/2.))*(math.log(r_x/sigma)))
    #      #V.append(epsilon*np.exp(-(r_x/sigma)**p))
    #      #V.append(epsilon*np.exp(-((r_x/sigma)**2)))
    #      #V.append(  epsilon * ( (sigma / r)**12))
    # ax1.plot(r,V)
    # print "Sigma:",sigma," p: ", p, " l_c: ", l_c

    ################
    f = filter(lambda x: 'hdf5' in x, os.listdir('.'))
    count = 0
    for filename in sorted(f):
        hpy = h5py.File(filename,'r')
        rbin = np.linspace(0,20,600)
        hist_total = np.zeros(rbin.shape[0]-1)
        for K in sorted(hpy.keys()):
            dist = hpy[K].value
            # Get length and Number of paricles
            N = hpy[K].attrs['N_atoms']
            L = hpy[K].attrs['L']
            num_frames = hpy[K].attrs["n_frames"]
            if 'ps2pe' in K:
                pass
            else:
                count+=1
                hist,rs = np.histogram(dist,bins=rbin)
                hist_total+=hist
        N = np.sum(hist_total)
        dr = rs[1]-rs[0]
        for i in range(1,hist_total.shape[0]):
            hist_total[i] = hist_total[i] / (count*num_frames*N*(4*math.pi*dr*rs[i]**2))
        scale = 1/hist_total[1:].max()
        ax2.plot(rs[:-1],Smooth(hist_total)*scale,'')
        ax2.set_ylabel('all')
        
        hpy.close

    sp4 = 7.5
    sp8 = 8
    sp12 = 8


if __name__ == "__main__":
    rmin=6
    rmax=20
    #full_radial_plots(ax1)
    #ax1.set_xlim(0,rmax*3)
    #ax1.set_ylim(0,3.5)
    #plt.show()
    gen_plots()
    ax1.legend(loc=1)
    #ax1.set_ylim(.01,5)

    ax1.set_yscale('log')
    ax1.set_ylim(.01,20)
    #ax2.set_ylim(.01,50)
    ax1.set_xlim(5,20)
    ax2.set_xlim(1.5,12)
    plt.show()


