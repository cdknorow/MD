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
def simpson_integration(f,x,initial=0):
    A = [initial]
    m = [x[0]]
    for i in range(0,x.shape[0]-2,2):
        A.append((x[i+2]-x[i])/6. * (f[i]+4*f[i+1]+f[i+2])+A[-1])
        m.append(x[i+2])
    return A, m


#autocorrelation function
#b is block number
#n_b is number of blocks
#tau_b is the number of steps in
# <A>_b = 1/tau_b SUM{tau=1}{tau_b}A(t) #sigma^2(<A>_b) = 1 / n_b SUM{b=1}{n_b}(<A>_b - <A>_run)^2
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

    # get the number of particles and box volume at each time step
    #N = int(os.getcwd().split('_')[14])
    N = 1000
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
    #del Volume[-1]
    #del breaks[-1]
    scale = 4.5
    #scale = 1.0
    for i in range(len(Volume)):
        print i , Volume[i]**(1/3.)/scale
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
    def statistical_iniefficiency():
        frames = [0]
        ax1.set_xlabel(r'$\tau_B^{1/2}$')
        ax1.set_ylabel('s')
        for i in frames:
            sp, x = autocorrelation(pressure[i])
            st, x = autocorrelation(temperature[i])
            p, x = autocorrelation(potential[i])
            ax1.plot(x,sp,'o',ms=10,label='pressure')
            ax1.plot(x,st,'^',ms=10,label='temperature')
            ax1.plot(x,p,'x',ms=10,label='potential')

    # Break down into averages
    #print "  Pressue:",P,"  Temperature:",T #print "  N * T / V = ", N*T/ V 
    ########################################################
    def PvsV():
        print 'pressure     potential        temperature         Volume' 
        for i in range(10):
            print ('%.4f        %.4f         %.4f        %.4f '%(pressure_avg[i], potential_avg[i], temperature_avg[i], Volume[i]))
        global counter
        ax1.set_xlabel(r'$\rho$')
        ax1.set_ylabel('P')
        for i in range(len(Volume)):
            NV = []
            for j in range(len(pressure[i])):
                NV.append(N/Volume[i])
            ax1.plot(NV, np.array(pressure[i])/max(pressure_avg),'-%c'%colors[counter])
            ax1.plot(NV, np.array(potential[i])/max(potential_avg),'-%c'%colors[counter])
        ax1.plot( N/np.array(Volume), np.array(pressure_avg)/max(pressure_avg),'x',label='pressure')
        ax1.plot( N/np.array(Volume), np.array(potential_avg)/max(potential_avg),'x',label='potential')



    ######################################################33
    def PVvsV_correlation():
        #fstart sphere
        #sigma=2.25
        #sigma = .5
        global counter
        TP = []
        NV = []
        Z = []
        Z_var = []
        print 'rho  Z      Z_var     Z_var/Z*100%'
        for i in range(len(Volume)):
            TP = []
            for j in range(len(temperature[i])):
                TP.append(pressure[i][j]*Volume[i]/(temperature[i][j]*N))
            NV.append(N*scale**3/(Volume[i]))
            Z.append(sum(TP)/len(TP))
            var = 0
            for i in range(len(TP)):
                var += (TP[i]- Z[-1])**2/len(TP)
            Z_var.append(var**0.5)
            print NV[-1],  Z[-1],  Z_var[-1]
        ax1.errorbar(NV,Z,yerr=Z_var,fmt='x%c'%colors[counter],label=title)
        ax1.set_xlabel(u'$\eta$')
        ax1.set_ylabel('PV/NkT')
        #ax1.set_ylim([1.0,50])
        #ax1.set_xlim([0.0001,.1])

    ######################################################33
    def P_correlation():
        global counter
        TP = []
        NV = []
        Z = []
        Z_var = []
        for i in range(len(Volume)):
            TP = []
            for j in range(len(temperature[i])):
                TP.append(pressure[i][j])
            NV.append(N*4*math.pi*(2.25)**3/(3*Volume[i]))
            Z.append(sum(TP)/len(TP))
            var = 0
            for i in range(len(TP)):
                var += (TP[i]- Z[-1])**2/len(TP)
            Z_var.append(var**0.5)
            print Z[-1], Z_var[-1], Z_var[-1]/Z[-1]
        ax1.errorbar(NV,Z,yerr=Z_var,fmt='x%c'%colors[counter],label=title)
        #plt.semilogx()
        #plt.semilogy()
        ax1.set_xlabel(u'$\eta$')
        ax1.set_ylabel('P')


    ########################################################
    #plot the change in free energy
    def FreeEnergy():
        global counter
        sigma=1.0
        gamma = [0]
        Z = [0]
        NV = [0]
        Z_potential = [0]
        T_avg = sum(temperature_avg)/len(temperature_avg)
        for i in range(len(Volume)):
            # rho01
            NV.append(N/(Volume[i]))
            # z - 1
            Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1)/NV[-1])
        del gamma[1]
        del NV[1]
        del Z[1]
        NV = np.array(NV)
        F_excess_Z, NV = simpson_integration(Z, NV,initial=0)
        print 'gamma_6      NV      Free Energy      Free Energy Potential      Z-1'
        for i in range(4):
            print ('%.4f            %.4f         %.4f         %.4f       %4f   '%(gamma[i], NV[i],
               F_excess_Z[i], F_excess_potential[i] , Z[i]))
        #print '   gamma_6       Volume    Pressure     Potential'
        print len(F_excess_Z)
        #plt.semilogx()
        #plt.semilogy()
        ax1.set_xlabel(u'$\gamma_{%i}$'%int(n*3))
        ax1.set_ylabel('Excess Free Energy of Liquid')
        #ax1.plot(gamma,Z_potential, '-x',label=title+' Z_potential')
        ax1.plot(NV,F_excess_Z, '-%c'%colors[counter],label=title+' Z-1')
        #ax1.plot(gamma,F_excess_potential, '-',label=title+' potential')
        #U6 potential (1992)


        ax1.set_ylim([0,30])
        ax1.set_xlim([0,2.7])
        #ax1.set_xlim([1.5,1.7])
        #ax1.set_ylim([18.3,18.4])
        #ax1.set_xlim([1.86,1.95])
        #ax1.set_ylim([1,8])
        #ax1.set_xlim([.3,1])
        #ax1.set_ylim([25,33])
        #ax1.set_xlim([2.2,2.7])

    if delta%2:
        #statistical_iniefficiency()
        PVvsV_correlation()
        #P_correlation()
    else:
        FreeEnergy()
    #PvsV()

 
    #statistical_iniefficiency()
    #PT_plot()
    
    #pot_increase()
    counter+=1



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
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    #ax2 = fig1.add_subplot(212)
    if len(arg) > 1:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = 'sp %i '%(int(f.split('_')[3])+1) + ' N '+f.split('_')[13]
            pvt(ax1, title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    elif len(arg) ==1:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = 'sp %i '%(int(f.split('_')[3])+1) + ' N '+f.split('_')[13]
            pvt(ax1, title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    else:
            title = 'title'
            pvt(ax1, title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
    plt.legend(loc=2)
    #import carnahan
    #reload(carnahan)
    #reload(peos)
    #carnahan.starling(ax1)
    #peos.modifiedstarling(ax1,c1=1,c2=2,c3=3)
    plt.show()

