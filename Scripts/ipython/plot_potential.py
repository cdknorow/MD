import os
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
import powerexpansioneos as peos
reload(peos)

if __name__ == '__main__':
    plt.close()
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    for file_potential in os.listdir('.'):
        print file_potential
        fid = open(file_potential,'r')
        V = []
        F = []
        r = []
        if file_potential[:3] == 'sp4':
            cut=11.5
            xmin=20
            sigma=8.68
            epsilon=1
            p=9.49
            l_c=1.28
            label = 'sp4'
        if file_potential[:3] == 'sp8':
            xmin=15
            cut=13.5
            sigma=10
            epsilon=1
            p=10
            l_c=1.2
            label = 'sp8'
        if file_potential[:4] == 'sp12':
            cut=15.5
            xmin=20
            sigma=13.6
            epsilon=1
            p=10.7
            l_c=1.45
            label = 'sp12'
        for line in fid.readlines():
            if line[0] != '#':
                x = float(line.split()[0])
                if  x < cut:
                    r.append(float(line.split()[0]))
                    V.append(float(line.split()[1]))
                    F.append(float(line.split()[2]))
        fid.close()
        ax1.plot(r,V,'-',label=label)
    ax1.set_yscale('log')
    #ax1.set_xscale('log')
    ax1.set_ylim((.01,50))
    ax1.set_xlim((6,17))
    #plt.legend(loc=1)
    ax1.set_xlabel('r')
    ax1.set_ylabel('U(r)')
    epsilon = 1.0 
    sigma = 18
    V=[]
    for r_x in r:
         V.append(-epsilon*5/18.*(30**(3/2.))*(math.log(r_x/sigma)))
    ax1.plot(r,V,label='-log(r/sigma)')
    plt.legend(loc=3)
    plt.show()
