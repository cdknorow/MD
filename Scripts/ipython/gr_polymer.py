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
matplotlib.rc('font', **{'size':'14'})
counter = 0
colors = ['r','g','b','k','y','c','m']
choose = [1,2,3,4]

def Smooth(M):
    #Find the max
    M_new = []
    for i in range(len(M)-1):
        M_new.append((M[i-1]+M[i]+M[i+1])/3.)
    M_new.append(M[-1])
    return np.array(M_new)


def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(311)
    ax2 = fig1.add_subplot(312)
    ax3 = fig1.add_subplot(313)
    return ax1, ax2, ax3

def radial_plots(ax1, ax2, ax3):
    hpy = h5py.File('polymer_stetching.hdf5','r')
    
    rbin = np.linspace(2,11,500)
    rbin2 = np.linspace(0.5,9,200)
    #rbin = np.linspace(5,12,500)
    #rbin2 = np.linspace(0.5,8,200)

    r = np.zeros(rbin.shape[0]-1)
    for i in range(1,rbin.shape[0]):
        r[i-1] = ((rbin[i]+rbin[i-1])/2)
    dr = r[1]-r[0]
    hist_total = np.zeros(r.shape)
    count = 0
    counter = 0

    for K in sorted(hpy.keys()):
        if counter in choose:
            dist = hpy[K].value
            # Get length and Number of paricles
            N = hpy[K].attrs['N_atoms']
            L = hpy[K].attrs['L']
            num_frames = hpy[K].attrs["n_frames"]
            if 'ps2pe' in K:
                hist,rs = np.histogram(dist,bins=rbin2)
            else:
                hist,rs = np.histogram(dist,bins=rbin)
                hist_total+=hist
            A = sum(hist)
            if 'ps2pe' in K:
                print counter, 'L =',L,'N =',N, "Dataset",K[:-5],
                print 'N_frames = ', num_frames
                ax3.plot(rs[:-1],Smooth(hist)/A,'')
                ax3.set_ylabel('ps2pe')
            if 'c2pe' in K:
                ax1.plot(r,Smooth(hist)/A,'')
                count+=1
            if 'c2ps' in K:
                ax1.plot(r,Smooth(hist)/A,'')
            if 'c2sp' in K:
                ax1.plot(r,Smooth(hist)/A,'')
                ax1.set_ylabel('c2sp')
        if 'ps2pe' in K:
            counter+=1
    N = np.sum(hist_total)
    for i in range(hist_total.shape[0]):
        hist_total[i] = hist_total[i] / (count*num_frames*N*(4*math.pi*dr*r[i]**2))
    ax2.plot(2*r,Smooth(hist_total),'')
    ax2.set_ylabel('all')
    
    hpy.close



if __name__ == "__main__":
    ax1,ax2,ax3 = get_plot()

    radial_plots(ax1,ax2,ax3)
    #ax2.set_xlim((5,11.5))
    #ax2.set_xlim((5,7.5))
    plt.legend(loc='upper right')
    plt.show()


