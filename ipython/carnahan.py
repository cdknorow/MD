import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
matplotlib.rc('font', **{'size':'20'})

def starling(ax1):
    Z_carnahan = []
    Z_exact6 = []
    Z_exact4 = []
    Z_exact10 = []
    Z_exact8 = []

    eta = np.linspace(0.0001,.45,num=1000)
    for x in eta:
        Z_carnahan.append((((1 + x + x**2 -x**3)/(1-x)**3) -1))
        Z_exact4.append(1 + 4*x + 10*x**2 + 18.365*x**3)
        Z_exact6.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5)
        Z_exact8.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
       533.34*x**6 + 68.54*x**7)
        Z_exact10.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
       533.34*x**6 + 68.54*x**7 + 85.8*x**8 + 105.8*x**9)

    ax1.plot(eta,Z_carnahan,'--',label='Carnahan-Starling')
    ax1.plot(eta,Z_exact4,label='4 terms')
    ax1.plot(eta,Z_exact6,label='6 terms')
    ax1.plot(eta,Z_exact8,label='8 terms')
    ax1.plot(eta,Z_exact10,label='10 terms')

if __name__ == '__main__':
    plt.close()
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    xlabel = r'$\eta$'
    ylabel = 'Z'

    Z_carnahan = []
    Z_exact6 = []
    Z_exact4 = []
    Z_exact10 = []
    Z_exact8 = []

    eta = np.linspace(0,.2,num=1000)
    #for x in eta:
    #    Z_carnahan.append((1 + x + x**2 -x**3)/(1-x)**3) 
    #    Z_exact4.append(1 + 4*x + 10*x**2 + 18.365*x**3)
    #    Z_exact6.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5)
    #    Z_exact8.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
    #   533.34*x**6 + 68.54*x**7)
    #    Z_exact10.append(1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
    #   533.34*x**6 + 68.54*x**7 + 85.8*x**8 + 105.8*x**9)

    for x in eta:
        Z_carnahan.append((((1 + x + x**2 -x**3)/(1-x)**3)-1)/(x*6/math.pi))
        Z_exact4.append((((1 + 4*x + 10*x**2 + 18.365*x**3))-1)/(x*6/math.pi))
        Z_exact6.append((((1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5))-1)/(x*6/math.pi))
        Z_exact8.append((((1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
       533.34*x**6 + 68.54*x**7))-1)/(x*6/math.pi))
        Z_exact10.append((((1 + 4*x + 10*x**2 + 18.365*x**3 + 28.224*x**4 + 39.82*x**5 +
       533.34*x**6 + 68.54*x**7 + 85.8*x**8 + 105.8*x**9))-1)/(x*6/math.pi))
    print Z_carnahan[1]
    plt.semilogy()
    plt.semilogx()
    ax1 = fig1.add_subplot(111)
    ax1.plot(eta,Z_carnahan,'--',label='Carnahan-Starling')
    ax1.plot(eta,Z_exact4,label='4 terms')
    ax1.plot(eta,Z_exact6,label='6 terms')
    ax1.plot(eta,Z_exact8,label='8 terms')
    ax1.plot(eta,Z_exact10,label='10 terms')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_ylim(1,15)
    ax1.set_xlim(.01,.6)

    ax1.legend(loc=2)
    plt.show()
