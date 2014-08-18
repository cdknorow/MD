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
choose = range(20)

def Smooth(M):
    #Find the max
    M_new = []
    for i in range(len(M)-1):
        M_new.append((M[i-1]+M[i]+M[i+1])/3.)
    M_new.append(M[-1])
    return np.array(M_new)


def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212)

    return ax1, ax2

def radial_plots(ax1, ax2):
    hpy = h5py.File('polymer_stetching.hdf5','r')
    for K in sorted(hpy.keys()):
        if counter in choose:
            dist = hpy[K].value
            # Get length and Number of paricles
            N = hpy[K].attrs['N_atoms']
            L = hpy[K].attrs['L']
            num_frames = hpy[K].attrs["n_frames"]
            hist,r = np.histogram(dist,bins=100,)
            r_new = np.zeros(r.shape[0]-1)
            for i in range(1,r.shape[0]):
                r_new[i-1] = ((r[i]+r[i-1])/2)
            r = r_new
            if 'ps2pe' in K:
                print counter, 'L =',L,'N =',N, "Dataset",K[:-5],
                print 'N_frames = ', num_frames
                ax1.plot(r,Smooth(hist),'')
            if 'c2pe' in K:
                ax2.plot(r,Smooth(hist),'')
            if 'c2ps' in K:
                ax2.plot(r,Smooth(hist),'')
    hpy.close


if __name__ == "__main__":
    ax1,ax2 = get_plot()

    radial_plots(ax1,ax2)
    ax2.set_xlim((5,11.5))
    #ax2.set_xlim((5,7.5))
    plt.legend(loc='upper right')
    plt.show()


