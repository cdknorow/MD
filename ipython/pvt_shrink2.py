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


# Integrate using simpsons rule
# int_a^b f(x) dx ~ (b-a)/6 * (f(a) + 4 * f((a+b)/2) + f(b))
def simpson_integration(f,x,initial=0):
    A = [initial]
    m = [x[0]]
    for i in range(0,len(x)-2,2):
        A.append((x[i+2]-x[i])/6. * (f[i]+4*f[i+1]+f[i+2])+A[-1])
        m.append(x[i+2])
    return np.array(A), np.array(m)

# Integrate using simpsons rule
# int_a^b f(x) dx ~ (b-a)/6 * (f(a) + 4 * f((a+b)/2) + f(b))
def simpson_integration_min(f,x,initial=0,int_min=.1):
    A = [initial]
    m = [int_min]
    for i in range(0,len(x)-2,2):
        if x[i+1]>int_min:
            A.append((x[i+2]-x[i])/6. * (f[i]+4*f[i+1]+f[i+2])+A[-1])
            m.append(x[i+2])
            #print i,x[i+2],A[-1]
    return np.array(A), np.array(m)

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


def pvt(title, xmin=0,xmax=1,delta=1):
    global counter
    global colors
    global n
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
    n = float(os.getcwd().split('/')[-1].split('_')[0][1:])
    print 'n is set to ',n
    #n = n/3.
    # get the number of particles and box volume at each time step
    try:
        fid = open('nsphere.txt','r')
        N = int(fid.readline().split()[0])
        fid.close()
    except:
        fid = open('hardsphere.hoomd','r')
        for i in range(7):
            fid.readline()
        line = fid.readline().split('=')[-1]
        N = int(line.split('e')[0])*10**(int(line.split('e')[-1]))
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
    #del Volume[-1]
    #del breaks[-1]
    fid.close()
    #get the pressure at each time step
    fid = open('mylog.log','r')
    M = fid.readlines()
    fid.close()
    pressure, pressure_avg = log(M[1:],M[0].split().index('pressure'),breaks,delta)
    temperature, temperature_avg = log(M[1:],M[0].split().index('temperature'),breaks,delta)
    potential, potential_avg = log(M[1:],M[0].split().index('potential_energy'),breaks,delta)
    sigma_p = []


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
   ########################################################
    def PT_plot():
        frames = [0]
        ax1.set_xlabel(r'$Time$')
        ax1.set_ylabel('pressure')
        for i in frames:
            x = range(len(pressure[i]))
            ax1.plot(x,pressure[i],'o',ms=10,label='pressure')
            ax1.plot(x,temperature[i],'^',ms=10,label='temperature')
            ax1.plot(x,potential[i],'^',ms=10,label='potential')

    # Break down into averages
    #print "  Pressue:",P,"  Temperature:",T #print "  N * T / V = ", N*T/ V 
    ########################################################
    def PvsV(pfit_coef = [-2.062305, 10.301981, -6.364244 ,10.012567], 
             pfit_no_coef = [-2.062305, 10.301669 ,-6.364247 ,10.012568],
             pfit_fe_coef = [-0.626243, 4.208299, -0.225634, 2.708868]):

        ax1 = get_plot()
        def pressure_exact(rho,a,c1,c2,c3):
            return (1+a)*rho**1+c1*rho**2+c2*rho**3+c3*rho**4

        def pressure_exact_rho(rho,a,c1,c2,c3):
            return rho+a*rho**2+c1*rho**3+c2*rho**4+c3*rho**5

        def pressure_derived(rho,a,c1,c2,c3):
            return rho+c1*rho**2+2*c2*rho**3+3*c3*rho**4

        #print 'pressure     potential        temperature         Volume' 
        #for i in range(10):
        #    print ('%.4f        %.4f         %.4f        %.4f '%(pressure_avg[i], potential_avg[i], temperature_avg[i], Volume[i]))
        global counter
        ax1.set_xlabel(r'$\rho$')
        ax1.set_ylabel('P')
        pfit = []
        pfit_no = []
        pfit_fe = []
        gamma = []
        for i in range(len(Volume)):
            NV = []
            for j in range(len(pressure[i])):
                NV.append(N/Volume[i])
            gamma.append(NV[-1])
            ax1.plot(NV, np.array(pressure[i]),'-%c'%colors[counter])
            pfit.append(pressure_exact(NV[-1] ,pfit_coef[0],pfit_coef[1],pfit_coef[2],pfit_coef[3]))
            pfit_no.append(pressure_exact(NV[-1], pfit_no_coef[0],pfit_no_coef[1],pfit_no_coef[2],pfit_no_coef[3]))
            pfit_fe.append(pressure_derived(NV[-1] , pfit_fe_coef[0], pfit_fe_coef[1], pfit_fe_coef[2], pfit_fe_coef[3]))

        for rho in gamma:
            a = pressure_exact(rho,pfit_coef[0],pfit_coef[1],pfit_coef[2],pfit_coef[3])
            c = pressure_exact(rho,pfit_no_coef[0],pfit_no_coef[1],pfit_no_coef[2],pfit_no_coef[3])
            b = pressure_derived(rho, pfit_fe_coef[0], pfit_fe_coef[1], pfit_fe_coef[2], pfit_fe_coef[3])
            print 'rho = ', rho,
            print ' Pfit-FE_fit = ', a - b,
            print ' % error= ', (a - b)/a
        #print c - b
        #print d - b
        ax1.plot( N/np.array(Volume), pfit,'x',label='pfit')
        #ax1.plot( N/np.array(Volume), pfit_no,'x',label='pfit no')
        ax1.plot( N/np.array(Volume), pfit_fe,'x',label='pfit fe')
        ax1.plot( N/np.array(Volume), np.array(pressure_avg),'x',label='pressure')
        #ax1.plot( N/np.array(Volume), np.array(potential_avg)/max(potential_avg),'x',label='potential')

    ########################################################
    def ZvsZpot_gamma(pfit_fe_coef=[], set_fe = True):
        if set_fe:
            if n == 6.25:
                pfit_fe_coef =[-0.123471,3.112000,3.320848,0.083388]
            if n == 7:
                pfit_fe_coef = [-0.051930 ,2.703564 ,3.246362,0.379012]
            if n == 8:
                pfit_fe_coef = [-0.112304,2.741619 ,2.576233, 0.984856]
            if n == 9:
                pfit_fe_coef = [-0.282307, 3.300624 ,1.358420, 1.786398]
            if n == 10:
                pfit_fe_coef = [-0.626243, 4.208299, -0.225634, 2.708868]
            if n == 11:
                pfit_fe_coef = [-0.764629, 5.146447 ,-1.897390, 3.667005]
            if n == 12:
                pfit_fe_coef = [-0.951805, 6.005817, -3.441567, 4.552962]
            if n == 13:
                pfit_fe_coef = [-2.138254, 9.811643, -7.846507,6.339112]
            if n == 14:
                pfit_fe_coef = [-1.857837,9.347776,-8.089134,6.768926]
        global counter
        #ax1 = get_plot()
        #ax1.set_xlabel(r'$\gamma_6$')
        #ax1.set_ylabel('Z')
        ax1 = 1
        gamma=[]
        NV = []
        Z = []
        Z_potential = []
        r_cut = float(os.getcwd().split('/')[-1].split('_')[-2])
        #r_cut=4.2
        print 'N is ',N
        Z_potential_no = []
        for i in range(len(Volume)):
            gamma.append(N/Volume[i]*temperature_avg[i]**(-3/n))
            # rho
            NV.append(N/(Volume[i]))
            # (Z - 1) / rho
            Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1))
            # (2 * Beta Phi_T / N ) / rho
            DE  = N*NV[-1]*2*math.pi/((n-3)*r_cut**(n-3))
            #print 'Delta DE', DE/potential_avg[i]
            Z_potential.append(n/3*((potential_avg[i]+DE)/(N*temperature_avg[i])))
            Z_potential_no.append(n/3*((potential_avg[i])/(N*temperature_avg[i])))

        fmin=0
        fmax=-1
        fit_Z = []
        fit_ZFE = []
        print "potential"
        fit_Z_potential = peos.viralexpansion_3(ax1, Z_potential[fmin:fmax], np.array(gamma[fmin:fmax]),'r',show=True)
        print "potential w/o correction"
        fit_Z_potential_no = peos.viralexpansion_3(ax1, Z_potential_no[fmin:fmax], np.array(gamma[fmin:fmax]),'r',show=True)
        print "Free Energy"
        print 'a=%f,c1=%f,c2=%f,c3=%f'%(pfit_fe_coef[0],pfit_fe_coef[1],pfit_fe_coef[2],pfit_fe_coef[3])
        #fit_Z_pressure = peos.viralexpansion_3(ax1, Z[fmin:fmax], np.array(gamma[fmin:fmax]),'r',show=True) 
        PvsV(fit_Z_potential,fit_Z_potential_no,pfit_fe_coef)



    ########################################################
    def Z_gamma():
        global counter
        ax1 = get_plot()
        ax1.set_xlabel(r'$\gamma$')
        ax1.set_ylabel('Z-1')
        r_cut = float(os.getcwd().split('/')[-1].split('_')[-2])
        NV = []
        Z_potential = []
        for i in range(len(Volume)):
            # rho
            NV.append(N/(Volume[i]))
            # Z - 1
            #Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1)/NV[-1])
            # 2 * Beta Phi_T / N 
            DE  = N*NV[-1]*2*math.pi/((n-3)*r_cut**(n-3))
            Z_potential.append(n/3*((potential_avg[i]+DE)/(N*temperature_avg[i]*NV[-1])))
            #Z_potential_no.append(n/3*((potential_avg[i])/(N*temperature_avg[i]*NV[-1])))

        #ax1.errorbar(NV,Z,yerr=Z_var,fmt='x%c'%colors[counter],label=title)
        #peos.modifiedstarling_fit_all(ax1,np.array(Z_potential),np.array(NV),colors[counter],show=True)
        peos.viralexpansion_3(ax1,np.array(Z_potential),np.array(NV),colors[counter],show=True)


        ########################################################
    def pot_increase():
        print title
        print 'power', n
        print N
        ax1.set_xlabel(r'$\gamma$')
        ax1.set_ylabel('potential/NKT')
        T_avg = sum(temperature_avg)/len(temperature_avg)
        print T_avg
        gamma=[]
        for i in range(len(Volume)):
            gamma.append(N/(2**0.5*Volume[i])*temperature_avg[i]**(-1./n))
        ax1.plot(gamma, (np.array(potential_avg)-potential_avg[0])/(N*T_avg),'x',label='phi-phi_0/Nkt')
        

    ######################################################33
    def PVvsV_correlation():
        global counter
        ax1 = get_plot()
        #ax1 = get_plot()
        TP = []
        NV = []
        Z = []
        Z_var = []
        for i in range(len(Volume)):
            TP = []
            for j in range(len(temperature[i])):
                TP.append(pressure[i][j]*Volume[i]/(temperature[i][j]*N))
            NV.append(N/(Volume[i]))
            Z.append(sum(TP)/len(TP))
            var = 0
            for i in range(len(TP)):
                var += (TP[i]- Z[-1])**2/len(TP)
            Z_var.append(var**2)
            #print Z[-1], Z_var[-1], Z_var[-1]/Z[-1]*100
        NV[0] = 0
        Z[0] = 1
        Z_var[0]=0

        ax1.errorbar(NV,Z,yerr=Z_var,fmt='x%c'%colors[counter],label=title)
        #peos.modifiedstarling_fit_all(ax1,np.array(Z),np.array(NV),colors[counter])
        peos.viralexpansion_3(ax1,np.array(Z),np.array(NV),colors[counter],show=True)
        #plt.semilogx()
        #plt.semilogy()
        ax1.set_xlabel(u'$\eta$')
        ax1.set_ylabel('PV/NkT')
        #ax1.set_ylim([1.0,50])
        #ax1.set_xlim([0.01,2.4])

    ########################################################
    ########################################################
    #plot the change in free energy
    def FreeEnergy():
        global counter
        global n
        # Fit to Z-1 from N/V = 0 to g
        if n==6.:
            c1=3.586511
            c2=6.1187979
            c1= 3.630785
            c2=5.919866
        if n==7.:
            c1= 3.139660
            c2 = 6.038275
            g=.1
        if n==8:
            c1 = 2.889490
            c2=5.810482
            g=.1
        if n==9:
            c1= 2.748766
            c2=5.428522
            g=.1
        if n ==10:
            c1= 2.098992
            c2=7.387928
            g=.2
        if n ==11.:
            c1= 2.434736
            c2=5.788083
            g=.1
        if n ==12.:
            a=2.549678
            c1= 3.993234
            c2=2.712819
            c3=3.469902
            g=.1
        if n ==13.:
            a=1.000000
            c1= 2.532850
            c2=3.437891
            c3=4.843509
            g=.1
        if n ==14.:
            a=1.000000
            c1= 2.482896
            c2=3.367217
            c3=4.901456
            g=.1
        if n == 6.25:
            a=1.000000
            c1= 3.525116
            c2=5.409153
            c3=1.469899
            g=.1
        try:
            init =a*g+.5*c1*g**2+c2/3.*g**3+c3/4.*g**4
        except:
            init =c1*g+c2/2.*g**2
        r_cut = float(os.getcwd().split('/')[-1].split('_')[-2])
        #r_cut=4.2
        NV = [0]
        Z = [0]
        Z_potential = [c1]
        Z_potential_no = [c1]
        gamma2 = []
        for i in range(temperature_avg):
            temperature_avg[i]=1.0
        for i in range(len(Volume)):
            # rho
            NV.append(N/(Volume[i]))
            # Z - 1
            Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1)/NV[-1])
            # 2 * Beta Phi_T / N 
            DE  = N*NV[-1]*2*math.pi/((n-3)*r_cut**(n-3))
            Z_potential.append(n/3*((potential_avg[i]+DE)/(N*temperature_avg[i]*NV[-1])))
            Z_potential_no.append(n/3*((potential_avg[i])/(N*temperature_avg[i]*NV[-1])))
            ##
            gamma2.append(NV[-1]/2**0.5*(1/temperature_avg[i])**0.5)
        Z_potential[0] = c1
        #peos.viralexpansion(ax1,np.array(Z_potential[0:5]),np.array(NV[0:5]),colors[counter],show=True)
        #peos.modifiedstarling_fit_all(ax1,np.array(Z_potential),np.array(NV),colors[counter])
        ax1 = get_plot()

        #Integrate To get the Free Energy 
        F_excess_Z, gamma = simpson_integration_min(Z, NV,
               initial=init,int_min=g)
        F_excess_potential, NV = simpson_integration_min(Z_potential, NV,
               initial=init,int_min=g)
        gamma = []

        for i in NV:
            gamma.append(i*(sum(temperature_avg)/len(temperature_avg))**(-3/6.))
        gamma = np.array(gamma)
        ax1.plot(gamma,F_excess_potential, '-',label=title+' FE Potential')

        #FIT OF Free Energy
        fmin= 16 # Min Cutoff of fit
        print 'fmin =',fmin, 'fit range', NV[fmin],'-',NV[-1]
        fitted = peos.viralexpansion_3(ax1, F_excess_potential[fmin:], np.array(gamma[fmin:]),'r',show=True)

        ## fitted function ##
        fit_FE = []
        for g in gamma:
            fit_FE.append(fitted[0]+fitted[1]*g+fitted[2]*g**2+fitted[3]*g**3)
        ax1.plot(gamma,fit_FE, '-x',label=title+' fit')


        #fcc calculated trvsst
        fcc_excess = np.zeros((gamma.shape[0]))
        energy_excess = np.zeros((gamma.shape[0]))
        if n == 6:
            a_fcc = 3.6131
            b_fcc = 3.6982
            g_L = 2.36
            g_fcc = 2.39
        if n == 6.25:
            a_fcc = 3.4247505
            b_fcc = 3.7851997
        if n == 7:
            a_fcc = 2.9754654
            b_fcc = 4.0185839
            g_L = 2.05
            g_fcc = 2.08
        if n == 8:
            a_fcc = 2.5402261
            b_fcc = 4.2753817
            g_L = 1.7
            g_fcc = 1.73
        if n == 9:
            a_fcc = 2.2083911
            b_fcc = 4.4821760
            g_L = 1.49
            g_fcc = 1.52
        if n == 10:
            a_fcc = 1.9388997
            b_fcc = 4.6487343
            g_L = 1.4
            g_fcc = 1.43
        if n == 11:
            a_fcc = 1.7118338
            b_fcc = 4.7823161
            g_L = 1.25
            g_fcc = 1.28
        if n == 12:
            a_fcc = 1.5164850
            b_fcc = 4.8884526
            g_L = 1.16
            g_fcc = 1.2
        if n == 13:
            a_fcc = 1.3461175
            b_fcc = 4.9714449
            g_L = 1.10
            g_fcc = 1.2
        if n == 14:
            a_fcc = 1.1964035
            b_fcc = 5.0346944
            g_L = 1.0
            g_fcc = 1.03

        for i,g in enumerate(gamma):
            fcc_excess[i] = a_fcc*g**(n/3.) + b_fcc + (n/2.)*math.log(g)
            energy_excess[i] = a_fcc*g**(n/3.) + 3
        ax1.plot(gamma, fcc_excess, '-',label='fcc phonon')
        ax1.set_xlabel(u'$\gamma_{%i}$'%int(n))
        ax1.set_ylabel('Excess Free Energy of Liquid')
        coex_fcc  =  a_fcc*g_fcc**(n/3.) + b_fcc + (n/2.)*math.log(g_fcc)
        coex_liq  = fitted[0]+fitted[1]*g_L+fitted[2]*g_L**2+fitted[3]*g_L**3
        print 'Df = ', coex_fcc - coex_liq
    ########################################################
    ########################################################
    #plot the change in free energy
    def FreeEnergy_verbose():
        global counter
        global n
        # U6
        if n==6.:
            c1=3.586511
            c2=6.1187979
            c1= 3.630785 
            c2=5.919866
        if n==7.:
            c1= 3.139660 
            c2 = 6.038275
            g=.1
        if n==8:
            c1 = 2.889490 
            c2=5.810482
            g=.1
        if n==9:
            c1= 2.748766 
            c2=5.428522
            g=.1
        if n ==10:
            c1= 2.098992 
            c2=7.387928
            g=.2
        if n ==11.:
            c1= 2.434736 
            c2=5.788083
            g=.1
        if n ==12.:
            a=2.549678
            c1= 3.993234
            c2=2.712819
            c3=3.469902
            g=.1
        if n ==13.:
            a=1.000000
            c1= 2.532850
            c2=3.437891
            c3=4.843509
            g=.1
        if n ==14.:
            a=1.000000
            c1= 2.482896
            c2=3.367217
            c3=4.901456
            g=.1
        if n == 6.25:
            a=1.000000
            c1= 3.525116
            c2=5.409153
            c3=1.469899
            g=.1
        try:
            init =a*g+.5*c1*g**2+c2/3.*g**3+c3/4.*g**4
        except:
            init =c1*g+c2/2.*g**2
        r_cut = float(os.getcwd().split('/')[-1].split('_')[-2])
        #r_cut=4.2
        NV = [0]
        Z = [0]
        Z_potential = [c1]
        Z_potential_no = [c1]
        gamma2 = []
        P_T = []
        P_T_no = []
        #fid = open('../pv0t%s.csv'%title,'w')
        #fid2 = open('../pvt%s.csv'%title,'w')
        #fid.write('gamma,gamma2,PV_0/NKT,Z_potential_V0,Z_potential_V0(w/ tail correction)\n')
        #fid2.write('gamma,gamma2,PV/NKT,Z_potential,Z_potential w/tail correction)\n')
        for i in range(temperature_avg):
            temperature_avg[i]=1.0
        for i in range(len(Volume)):
            # rho
            NV.append(N/(Volume[i]))
            # Z - 1
            Z.append((pressure_avg[i]*Volume[i]/(N*temperature_avg[i])-1)/NV[-1])
            # 2 * Beta Phi_T / N 
            DE  = N*NV[-1]*2*math.pi/((n-3)*r_cut**(n-3))
            Z_potential.append(n/3*((potential_avg[i]+DE)/(N*temperature_avg[i]*NV[-1])))
            Z_potential_no.append(n/3*((potential_avg[i])/(N*temperature_avg[i]*NV[-1])))
            ##
            gamma2.append(NV[-1]/2**0.5*(1/temperature_avg[i])**0.5)
            #fid.write('%.4f,%.4f,%.4f,%.4f,%.4f\n'%(NV[-1]/temperature_avg[i],NV[-1]/2**0.5*(1/temperature_avg[i])**.25, 
            #                    pressure_avg[i]/(2**0.5*temperature_avg[i]), 
            #                     (Z_potential_no[-1]*NV[-1]+1)*NV[-1]/2**0.5,
            #                     (Z_potential[-1]*NV[-1]+1)*NV[-1]/2**0.5))

            #fid2.write('%.4f,%.4f,%.4f,%.4f,%.4f\n'%(NV[-1]/temperature_avg[i],NV[-1]/2**0.5*(1/temperature_avg[i])**.25, 
            #                    pressure_avg[i]/(temperature_avg[i]*NV[-1]),
            #                    (Z_potential_no[-1]*NV[-1]+1), 
            #                     (Z_potential[-1]*NV[-1]+1)))
            #P_T.append( ( (potential_avg[i]+DE) / (N*temperature_avg[i])) - 7.227*gamma2[-1]**2)
            #P_T_no.append( ( (potential_avg[i]) / (N*temperature_avg[i])) - 7.227*gamma2[-1]**2)
        #fid.close()
        #fid2.close()
        Z_potential[0] = c1
        #peos.viralexpansion(ax1,np.array(Z_potential[0:5]),np.array(NV[0:5]),colors[counter],show=True)
        #peos.modifiedstarling_fit_all(ax1,np.array(Z_potential),np.array(NV),colors[counter])
        #plt.show()
        ax1 = get_plot()

        #Integrate To get the Free Energy 
        F_excess_Z, gamma = simpson_integration_min(Z, NV,
               initial=init,int_min=g)
        F_excess_potential, NV = simpson_integration_min(Z_potential, NV,
               initial=init,int_min=g)
        gamma = []
        for i in NV:
            gamma.append(i*(sum(temperature_avg)/len(temperature_avg))**(-3/6.))
        gamma = np.array(gamma)
        fid=open('../FE%s.csv'%title,'w')
        fid.write('gamma,rho,FE_Z,FE_potential\n')
        for i in range(gamma.shape[0]):
            fid.write('%.4f,%.4f,%.4f,%.4f\n'%(gamma[i], NV[i], F_excess_Z[i], F_excess_potential[i]))
        fid.close()
        ax1.plot(NV,F_excess_Z, '-%c'%colors[counter],label=title+' FE Z-1')
        ax1.plot(gamma,F_excess_potential, '-',label=title+' FE Potential')

        #FIT OF Free Energy
        fmin= 16 # Min Cutoff of fit
        print 'fmin =',fmin, 'fit range', NV[fmin],'-',NV[-1]
        fitted = peos.viralexpansion_3(ax1, F_excess_potential[fmin:], np.array(gamma[fmin:]),'r',show=True)

        ## fitted function ##
        fit_FE = []
        fit_ZFE = []
        for g in gamma:
            #fit_ZFE.append(c1*g+c2/2.*g**2)
            fit_ZFE.append(a*g+.5*c1*g**2+c2/3.*g**3+c3/4.*g**4)
            fit_FE.append(fitted[0]+fitted[1]*g+fitted[2]*g**2+fitted[3]*g**3)
            #fit_FE.append(fitted[0]+fitted[1]*g+fitted[2]*g**2)
        ax1.plot(gamma,fit_FE, '-x',label=title+' fit')
        #ax1.plot(gamma,fit_ZFE, '-x',label=title+' Z fit')


        #fcc calculated trvsst
        fcc_excess = np.zeros((gamma.shape[0]))
        energy_excess = np.zeros((gamma.shape[0]))
        if n == 6:
            a_fcc = 3.6131
            b_fcc = 3.6982
        if n == 6.25:
            a_fcc = 3.4247505
            b_fcc = 3.7851997
        if n == 7:
            a_fcc = 2.9754654
            b_fcc = 4.0185839
        if n == 8:
            a_fcc = 2.5402261
            b_fcc = 4.2753817
        if n == 9:
            a_fcc = 2.2083911
            b_fcc = 4.4821760
        if n == 10:
            a_fcc = 1.9388997
            b_fcc = 4.6487343
        if n == 11:
            a_fcc = 1.7118338
            b_fcc = 4.7823161
        if n == 12:
            a_fcc = 1.5164850
            b_fcc = 4.8884526
        if n == 13:
            a_fcc = 1.3461175
            b_fcc = 4.9714449
        if n == 14:
            a_fcc = 1.1964035
            b_fcc = 5.0346944

        for i,g in enumerate(gamma):
            fcc_excess[i] = a_fcc*g**(n/3.) + b_fcc + (n/2.)*math.log(g)
            energy_excess[i] = a_fcc*g**(n/3.) + 3

        #bcc
        #bcc_excess = np.zeros((gamma.shape[0]))
        #a_bcc = 3.6303
        #b_bcc = 3.5647
        #for i,g in enumerate(gamma):
            #bcc_excess[i] = a_bcc*g**(n/3.) + b_bcc + (n/2.)*math.log(g)
        if counter in [0,1,2,3,4]:
            #print  'fcc bcc coexistance at gamma = %f'%((b_fcc - b_bcc)/(a_bcc-a_fcc))**0.5 
            #bcc and fcc excess free energies
            #ax1.plot(gamma, bcc_excess, '-',label='bcc exact')
            #ax1.plot(gamma, fcc_excess, '-',label='fcc exact')
            #bcc and fcc excess free energies - Simulation data for excess free energy of liquid
            #ax1.plot(gamma, bcc_excess-F_excess_potential, '-x',label='bcc_t - L_s')
            #ax1.plot(gamma, fcc_excess-F_excess_potential, '-x',label='fcc_t - L_s')
            #ax1.set_xlim([1,2.7])
            #ax1.set_ylim([-.1,1])c1= 2.753180 c2=6.438427
            pass
        paper_gamma = np.array([0.3162,0.6324,0.9487,1.2649,1.5811,1.8974,2.2136,2.3717])
        paper_liquid_FE = np.array([1.4609,3.5364,6.2504,9.6153,13.642,18.3382,23.7067,26.6442])
        RHNC_FE = np.array([1.4644,3.5407,6.2541,9.6180,13.6414,18.3307,23.6907,26.6233])
        sim_liquid_fit = []
        print 'gamma  fitted'
        #for i in paper_gamma:
            #sim_liquid_fit.append(fitted[0]+fitted[1]*i+fitted[2]*i**2+fitted[3]*i**3)
            #print i,sim_liquid_fit[-1]
        #sim_liquid_fit = np.array(sim_liquid_fit)
        bcc_excess = np.zeros(paper_gamma.shape)
        fcc_excess = np.zeros(paper_gamma.shape)
        for i,g in enumerate(paper_gamma):
            bcc_excess[i] = a_bcc*g**2 + b_bcc + (3)*math.log(g)
            fcc_excess[i] = a_fcc*g**2 + b_fcc + (3)*math.log(g)
        if counter==0:
            #pass
            #excess free energy for the liquid calculated from the paper
            #ax1.plot(paper_gamma,paper_liquid_FE, '-o',label='N=500 (1992)')
            #ax1.plot(paper_gamma,RHNC_FE, '-o',label=title+' RHNC')
            #ax1.plot(paper_gamma,sim_liquid_fit,'-^',label='Fitted')
            #ax1.plot(paper_gamma,paper_liquid_FE- sim_liquid_fit, '-o',label='N=500 (1992)')
            #ax1.plot(paper_gamma,RHNC_FE - sim_liquid_fit, '-o',label=title+' RHNC')
            #ax1.plot(paper_gamma[:-1],paper_liquid_FE[:-1]- F_excess_potential, '-o',label='N=500 (1992)')
            #ax1.plot(paper_gamma[:-1],RHNC_FE[:-1] - F_excess_potential, '-o',label=title+' RHNC')
            #ax1.set_xlim([1.1,2.5])
            #ax1.set_ylim([-.05,.05])
            #fcc excess free energy - liquid calculated from the paper
            #ax1.plot(paper_gamma,fcc_excess-paper_liquid_FE, '->',label='fcc - liquid(1992)')
            #ax1.plot(paper_gamma,fcc_excess-RHNC_FE, '-o',label='fcc_t - RHNC')
            #bcc excess free energy - liquid calculated from the paper
            #ax1.plot(paper_gamma,bcc_excess-paper_liquid_FE, '->',label='bcc_t - L_(1992)')
            #ax1.plot(paper_gamma,bcc_excess-RHNC_FE, '-o',label='bcc_t - RHNC')
            #ax1.set_xlim([1.5,2.5])
            #ax1.set_ylim([-.05,.05])
            pass
        paper_gamma = np.array([2.2136,2.2768,2.3401,2.4033])
        paper_fcc = np.array([23.7794,24.8952,26.0359,27.2020])
        paper_bcc = np.array([23.7668,24.8864,26.0321,27.2039])
        bcc_excess = np.zeros(paper_gamma.shape)
        fcc_excess = np.zeros(paper_gamma.shape)
        for i,g in enumerate(paper_gamma):
            bcc_excess[i] = a_bcc*g**2 + b_bcc + (3)*math.log(g)
            fcc_excess[i] = a_fcc*g**2 + b_fcc + (3)*math.log(g)
        if counter==0:
            #pass
            #excess free energy for the liquid calculated from the paper
            #ax1.plot(paper_gamma,bcc_excess-paper_bcc, '-o',label='N500 bcc - bcc_t' )
            #ax1.plot(paper_gamma,fcc_excess-paper_fcc, '-o',label='N500 fcc - fcc_t')
            #ax1.set_xlim([2,2.5])            #ax1.set_ylim([-.1,.1])
            pass
        if counter in [0]:
            pass
            #potential energy over static lattice (Phi_t - Phi_0)/NkT
            #gamma = gamma_6/(sqr(2))
            paper_gamma2 = np.array([.1,.25,.5,1,1.5])
            paper_fluid2 = np.array([.247,.577,1.036,1.739,2.195])
            #ax1.plot(paper_gamma2,paper_fluid2, '-o',label='N256 phi_T-Phi_e' )
            #ax1.plot(gamma2,P_T, '-o',label='Phi_T-Phi_e')
            #ax1.plot(gamma2,P_T_no, '-o',label='Phi_Tno-Phi_e')
            #ax1.set_xlabel(u'$\gamma_{%i}/(2^{.5})$'%int(n*3))
            #ax1.set_ylabel(u'$(\Phi_T-\Phi_0)/NkT$')
        ax1.set_xlabel(u'$\gamma_{%i}$'%int(n))
        ax1.set_ylabel('Excess Free Energy of Liquid')



    if delta%2:
        #statistical_iniefficiency()
        ZvsZpot_gamma()
        #Z_gamma()
        #PVvsV_correlation()
        #PvsV()
    else:
        #ZvsZpot_gamma()
        FreeEnergy()
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
    if len(arg) > 1:
        ax1 =get_plot()
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0]
            pvt(title = title,ax1 = ax1,xmin = int(config['-x']), xmax = int(config['-r']),
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
            pvt(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    plt.legend(loc=2)
    #import carnahan
    #reload(carnahan)
    #reload(peos)
    #carnahan.starling(ax1)
    #peos.modifiedstarling(ax1,c1=1,c2=2,c3=3)
    plt.show()


