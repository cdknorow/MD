import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate as integrate 
import math
import powerexpansioneos as peos
reload(peos)
#set the fonts
matplotlib.rc('font', **{'size':'20'})
counter = 0
colors = ['r','g','b','k','y','c','m']


def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def WHAM(title,xmin=0,xmax=1,delta=1):
    global counter
    global colors
    def log(M, index):
        analyze = []
        for i in range(len(M)):
            analyze.append(float(M[i].split()[index]))
        return analyze
    def log_distance(M, L):
        def boundary(x1,x2,L):
            if x1-x2 > L/2:
                return ((x1-x2)-L)**2
            if x1-x2 < -L/2:
                return ((x1-x2)+L)**2
            return (x1-x2)**2
        analyze = []
        for i in range(len(M)):
            line = M[i].split()
            x1 = float(line[2])
            y1 = float(line[3])
            z1 = float(line[4])
            x2 = float(line[5])
            y2 = float(line[6])
            z2 = float(line[7])
            r2 = boundary(x1,x2,L[0])+boundary(y1,y2,L[1])+boundary(z1,z2,L[2])
            analyze.append(r2**.5)
        return analyze
    fid = open('box_length.txt','r')
    M = fid.readline()
    L = [float(M.split()[0]), float(M.split()[1]), float(M.split()[2])]
    fid.close()
    ###################################################################
    fid = open('mylog.log','r')
    M = fid.readlines()
    fid.close()
    try:
        potential = log(M[1:],M[0].split().index('pair_lj_energy'))
    except:
        potential = log(M[1:],M[0].split().index('pair_table_energy'))
    fid = open('scry_potential.dat','r')
    M = fid.readlines()
    fid.close()
    r0 = log(M[1:],M[0].split().index('d'))


    r_set = [r0[0]]
    Ni = []
    count = 1
    for i in range(1,len(r0)):
        if r0[i] != r0[i-1]:
            r_set.append(r0[i])
            Ni.append(count)
            count = 0
        count+=1
    Ni.append(Ni[-1])
    r0 = r_set
    ######################################################
    distance = log_distance(M[1:],L)
    #distacnce_sim = log(M[1:],M[0].split().index('r'))
    potential_harmonic = []
    ax1 = get_plot()
    print r0
    plt.plot(distance,potential,'x')
    ax1 = get_plot()
    count = 0
    h =[]
    print Ni
    for i in range(len(distance)):
        if count % Ni[0]==0:
            hist, bin_edges = np.histogram(h,np.linspace(.8,2.5,800))
            x = []
            for i in range(bin_edges.shape[0]-1):
                x.append((bin_edges[i]+bin_edges[i+1])/2)
            ax1.plot(x,hist)
            h = []
        h.append(distance[i])
        count+=1
    plt.show()


    ax1 = get_plot()
    bins = np.linspace(.8,2.5,800)
    hist, bin_edges = np.histogram(distance,bins)
    x = []
    for i in range(bin_edges.shape[0]-1):
        x.append((bin_edges[i]+bin_edges[i+1])/2)
    ax1.plot(x,hist)

    plt.show()

    #Potential Bias U(x)_i
    print 'k = 50'
    def potential_bias(r, r0, k = 50):
        return k/2*(r-r0)**2

    # input
    #
    # Ns - Number of simulations
    # Ni - number of data points in the ith simulation
    # hist - histogram over x for all simulations
    # F - free energy shift for each simulation
    # x - x values to compute pmf for
    # r0 - biased potential r0 values
    def Px(Ni,  hist, x, r0, F, T=1.0):
        P = []
        for i in range(len(x)):
            B = 0
            for index in range(len(r0)):
                B += Ni[index] * math.exp((F[index] - potential_bias(x[i],r0[index]))/T)
            P.append(hist[i]/B)
        return P

    # input
    # Px - best estimate of ubiased probability distribution
    # x - x values to compute pmf for
    # r0 - biased potential r0 values
    def Fi(Px, x, r0, T=1.0):
        F = []
        for i in range(len(r0)):
            f = 0
            for j in range(len(x)):
                f += Px[i] * math.exp(-potential_bias(x[j],r0[i])/T)
            F.append(-T*np.log(f))
        return F

    #Guess the first F values
    F = []
    import random
    for i in range(len(r0)):
        F.append(0)
    ax1 = get_plot()
    #Run Iteration Until Self-Consistent
    for i in range(50):
        P = Px(Ni, hist, x, r0, F)
        F_new = Fi(P, x, r0)
        #for i in range(len(F_new)):
        #    print F_new[i],F[i]
        #print 'next iteration'
        F = F_new

    #ax1.plot(x[:200],P[:200])
    #plt.show()


    #ax1.plot(distance_file,potential_harmonic,'x')



    #get sigma^2



if __name__ == "__main__":
    import getopt            ## fitted function ##
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    plt.close()
    config = {'-x':'-1','-r':'-1','-i':0}
    options, arg =  getopt.getopt(sys.argv[1:],'x:r:D:i:')
    config.update( dict(options))
    if len(arg) > 1:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[-1]
            WHAM(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    else:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] 
            WHAM(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    plt.legend(loc=3)
    #import carnahan
    #reload(carnahan)
    #reload(peos)
    #carnahan.starling(ax1)
    #peos.modifiedstarling(ax1,c1=1,c2=2,c3=3)
    plt.show()

