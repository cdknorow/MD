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
# Integrate using simpsons rule
# int_a^b f(x) dx ~ (b-a)/6 * (f(a) + 4 * f((a+b)/2) + f(b))
def simpson_integration(f,x):
    A = []
    m = []
    for i in range(len(x)-2):
        A.append((x[i]-x[i+2])/6 (f[i]+4*f[i+1]+f[i+2]))
        m.append(x[i+1])
    return A, m




#autocorrelation function
#b is block number
#n_b is number of blocks
#tau_b is the number of steps in
# <A>_b = 1/tau_b SUM{tau=1}{tau_b}A(t)
#sigma^2(<A>_b) = 1 / n_b SUM{b=1}{n_b}(<A>_b - <A>_run)^2
def autocorrelation(A):
    tau = 50
    A_avg = sum(A)/len(A)
    #statistical inefficienty
    s = []
    sigma = 0
    #sigma(A)^2
    for t in range(len(A)):
        sigma += (A[t]-A_avg)**2
    sigma_A = sigma/len(A)
    #print A[10], A_avg, sigma_A
    #sigma(A)b^2
    for tau_b in range(1,tau,5):
        A_b = []
        n_b = len(A)/tau_b
        #print tau_b,n_b
        sigma = 0
        for b in range(n_b):
            Asum = 0
            for i in range(tau_b*b,tau_b*b+tau_b):
                Asum += A[i]
            sigma += (Asum/tau_b - A_avg)**2
        # simga^2(<A>_b)
        sigma_Ab =  (sigma/n_b)
        #print sigma_A, sigma_Ab
        s.append(tau_b * sigma_Ab/sigma_A)
    return s,np.array(range(1,tau,5))**0.5


def pvt(ax1, title,xmin=0,xmax=1,delta=1):
    global counter
    global colors
    pvt_avg = [[],[]]
    def log(M, index, breaks, delta):
        analyze = [[] for i in range(len(breaks))]
        analyze_avg = []
        count = 0
        for k in range(len(breaks)):
            count_xmin = 0
            for i in range(breaks[k][0],breaks[k][1]):
                if count_xmin > delta:
                    analyze[k].append(float(M[i].split()[index]))
                count_xmin += 1
            analyze_avg.append(sum(analyze[k])/(len(analyze[k])))
        return analyze, analyze_avg
    #get the exponent
    try:
        n = float(title.split('U')[-1])/3   
    except:
        n = float(title[1:3])/3   
    # get the number of particles and box volume at each time step
    fid = open('nsphere.txt','r')
    N = int(fid.readline().split()[0])
    fid.close()
    fid = open('mylog.log','r')
    fid.readline()
    V = []
    breaks = []
    Volume = []

    for line in fid.readlines():
        V.append(float(line.split()[4]))
    i = 0
    while i < len(V)-2:
        if V[i] == V[i+2]:
            dt = V.count(V[i])
            breaks.append([i,i+dt-1])
            Volume.append(V[i])
            i+=dt
        else:
            i+=1
    breaks = breaks[xmin:xmax]
    Volume = Volume[xmin:xmax]
    del Volume[-1]
    del breaks[-1]
    fid.close()
    #get the pressure at each time step
    fid = open('pressure.log','r')
    P = fid.readlines()
    fid.close()
    pressure, pressure_avg = log(P[1:],1,breaks,delta)
    sigma_p = []

    #get the temperature at each time step
    fid = open('mylog.log','r')
    T = fid.readlines()
    fid.close()
    temperature, temperature_avg = log(T[1:],1,breaks,delta)

    #get the potential energy at each time step
    fid = open('mylog.log','r')
    P = fid.readlines()
    fid.close()
    potential, potential_avg = log(P[1:],2,breaks,delta)
    #get sigma^2


    ########################################################
    ########################################################
    #plot the change in free energy
    def Get_Z():
        global counter
        sigma=1.0
        gamma = []
        NV = []
        Z = []
        Z_potential = []
        T_avg = sum(temperature_avg)/len(temperature_avg)
        for i in range(len(Volume)):
            # gamma 
            gamma.append(N/(Volume[i])*temperature_avg[i]**(-1./n))
            # rho01
            NV.append(N/(Volume[i]))
            # z - 1
            Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1)/NV[-1])
            # 2 * Beta Phi / N 
            Z_potential.append(n*potential_avg[i]/(N*temperature_avg[i]*NV[-1]))
        return [gamma,NV,Z]
    A = Get_Z()
    return A





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
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    #ax2 = fig1.add_subplot(212)
    gvz= []
    if len(arg) > 1:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0]
            gvz.append(pvt(ax1, title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i'])))
            os.chdir('../')
    plt.legend(loc=2)
    n = 6/3.0
    #get the parts of gvz that we want from each segment
    cut=-1
    gamma=[0]
    Z=[0]
    NV=[0]
    gamma += gvz[0][0][15:cut]
    NV += gvz[0][1][15:cut]
    Z += gvz[0][2][15:cut]
    cut = 5
    #print len(gamma)

    gamma += gvz[1][0][cut:]
    NV += gvz[1][1][cut:]
    Z += gvz[1][2][cut:]
    for i in range(len(gamma)):
        print i, gamma[i]

    F_excess_Z = simpson_integration(Z,NV)
    gamma_new = []
    for i in range(len(gamma)):
        if i%2:
            gamma_new.append(gamma[i])
    gamma = np.array(gamma_new)
    #print 'gamma_6  NV Free Energy   Z-1'
    #for i in range(len(gamma)-1):
    #    print ('%.4f            %.4f         %.4f         %.4f     '%(gamma[i], NV[i],
    #     F_excess_Z[i], Z[i]))
    #print '   gamma_6       Volume    Pressure     Potential'
    #plt.semilogx()
    ax1.set_xlabel(u'$\gamma_{%i}$'%int(n*3))
    ax1.set_ylabel('Excess Free Energy of Liquid')
    #ax1.plot(gamma,Z_potential, '-x',label=title+' Z_potential')
    ax1.plot(gamma,F_excess_Z, '-x',label=title+' Z-1')
    #ax1.plot(gamma,F_excess_potential, '-',label=title+' potential')
    #U6 potential (1992)
    real_data_gamma = [0.3162,0.6324,0.9487,1.2649,1.5811,1.8974,2.2136,2.3717]
    real_data_FE = [1.4609,3.5364,6.2504,9.6153,13.642,18.3382,23.706,26.644]
    ax1.plot(real_data_gamma,real_data_FE, 'o',label=title+' Paper Values')
    #bcc calculated trvsst
    gamma = np.linspace(1,4)
    bcc_excess = np.zeros((gamma.shape[0]))
    fcc_excess = np.zeros((gamma.shape[0]))
    a_bcc = 3.6303
    b_bcc = 3.5647
    a_fcc = 3.6131
    b_fcc = 3.6982
    for i,g in enumerate(gamma):
        bcc_excess[i] = a_bcc*g**2 + b_bcc + 3*math.log(g)
        fcc_excess[i] = a_fcc*g**2 + b_fcc + 3*math.log(g)
    print  'coexistance at gamma = %f'%((b_fcc - b_bcc)/(a_bcc-a_fcc))**0.5 
    ax1.plot(gamma, bcc_excess, '-',label='bcc excess')
    ax1.plot(gamma, fcc_excess, '-',label='fcc excess')
    ax1.set_ylim([0,5])
    ax1.set_xlim([.001,.7])
    plt.show()
   