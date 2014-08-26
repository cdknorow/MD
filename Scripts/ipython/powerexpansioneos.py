import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
matplotlib.rc('font', **{'size':'20'})

def modifiedstarling(ax1,c1=1,c2=2,c3=3):
    Z = []

    eta = np.linspace(0.0001,.45,num=1000)
    for x in eta:
        Z.append(c1*(1/(1-x)) + c2*x/(1-x)**2+c3*x**2/(1-x)**3)

    ax1.plot(eta,Z,label='C modified')

def modifiedstarling_fit_all(ax1, Z, eta,color,show=True):
    def coeff(x,a,c2,c3):
            return 1.0*(1/(1-a*x)) + c2*a*x/(1-a*x)**2+c3*(a*x)**2/(1-a*x)**3
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'alpha=%f c1= 1.0 c2=%f c3=%f'%(popt[0], popt[1],popt[2])
    if show:
        ax1.plot(eta,coeff(eta,popt[0],popt[1],popt[2]),
                    color )
        ax1.plot(eta,Z, 'x' )
        print 'difference Zactual-Zfit'
        for i in range(len(eta)):
            print eta[i], Z[i] - (1.0*(1/(1-popt[0]*eta[i])) + 
                    popt[1]*popt[0]*eta[i]/(1-popt[0]*eta[i])**2+
                    popt[2]*(popt[0]*eta[i])**2/(1-popt[0]*eta[i])**3)

def modifiedstarling_fit_alpha(ax1, Z, eta,color):
    def coeff(x,a):
            return 1.0*(1/(1-a*x)) + 2*a*x/(1-a*x)**2+3*(a*x)**2/(1-a*x)**3
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'alpha=%f c1= 1.0 c2=2 c3=3'%(popt[0])
    ax1.plot(eta,coeff(eta,popt[0]), color )

def modifiedstarling_fit_c(ax1, Z, eta,color):
    def coeff(x,c1,c2,c3):
            a=1.0
            return 1*(1/(1-a*x)) + c2*a*x/(1-a*x)**2+c3*(a*x)**2/(1-a*x)**3
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'alpha=1.0 c1= %f c2=%f c3=%f'%(popt[0], popt[1],popt[2])
    ax1.plot(eta,coeff(eta,popt[0],popt[1],popt[2]),
                colorw)

def viralexpansion_fit(ax1, Z, eta, color, show = False):
    def coeff(x,c1,c2):
            a=1.0
            return 1 + c1*x + c2*x**2
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'alpha=1.0 c1= %f c2=%f'%(popt[0], popt[1])
    g=.2
    c1=popt[0]
    c2=popt[1]
    g=.2
    print 'FE g=%f'%g,c1*g+c2/2.*g**2
    g=.3162
    print 'FE g=%f'%g,c1*g+c2/2.*g**2
    if show:
        ax1.plot(np.arange(eta[0],eta[-1],.01),coeff(np.arange(eta[0],eta[-1],.01),popt[0],popt[1]),
                color )
        ax1.plot(eta,Z, 'x' )
    print 'difference Zactual-Zfit'
    for i in range(len(eta)):
        print eta[i], Z[i] - (1+popt[0]*eta[i]+popt[1]*eta[i]**2)

def viralexpansion(ax1, Z, eta, color, show=False):
    def coeff(x,a,c1,c2):
            return a + c1*x + c2*x**2 
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'a=%f, c1= %f, c2=%f'%(popt[0], popt[1],popt[2])
    a=popt[0]
    c1=popt[1]
    c2=popt[2]
    if show:
        ax1.plot(np.linspace(0,eta[-1],100),coeff(np.linspace(0,eta[-1],100),popt[0],popt[1],popt[2]),
                color )
        ax1.plot(eta,Z,'x')
    for i in range(len(eta)):
        dx =  (Z[i] - (popt[0]+popt[1]*eta[i]+popt[2]*eta[i]**2))/Z[i]
        if abs(dx*100) > .1:
            print eta[i],abs(dx*100)

    return a,c1,c2

def viralexpansion_3(ax1, Z, eta, color, show = False):
    def coeff(x,a,c1,c2,c3):
        return a + c1*x + c2*x**2  + c3*x**3
    popt, pcov = curve_fit(coeff, eta, Z)
    print 'a=%f, c1= %f, c2=%f, c3=%f'%(popt[0], popt[1],popt[2],popt[3])
    print 'range %f - %f'%(eta[0],eta[-1])
    if show:
        fit_x = np.linspace(eta[0],eta[-1],100)
        ax1.plot(fit_x,coeff(fit_x,popt[0],popt[1],popt[2],popt[3]),
                color )
        ax1.plot(eta,Z,'x')
    for i in range(len(eta)):
        dx = (Z[i] - (popt[0]+popt[1]*eta[i]+popt[2]*eta[i]**2+popt[3]*eta[i]**3))/Z[i]
        if abs(dx*100) > .001:
            print eta[i],abs(dx*100)

    return [popt[0],popt[1],popt[2],popt[3]]

def Inverse_Potential(ax1,x,U,show=False):
    def coeff(r, sigma, p, epsilon):
        return epsilon * ( sigma/r)**p

    popt, pcov = curve_fit(coeff, x, U, p0 = [2.5, 14, 1], maxfev=50000)
    print 'sigma= %f p=%f epsilon%f'%(popt[0], popt[1],popt[2])
    if show:
        xmin=1.3
        xmax=4
        ax1.plot(x,U,'xr')
        ax1.plot(np.linspace(xmin,xmax,50),coeff(np.linspace(xmin,xmax,50),popt[0],popt[1],popt[2]),'-g')
    return [popt[0], popt[1],popt[2]]

def Soft_Potential(ax1,x,U,show=False):
    l_c = 1
    epsilon = 1
    def coeff(r, sigma, p, k, b):
        c = math.log(1+2**0.5)
        A =  1 / (np.sinh((r/sigma))**p)
        #w = r / (sigma * l_c)
        #B = k * np.cos(math.pi*((1-b)*w / (b+(1-2*b)*w)))
        #B = k*np.exp(-(r-sigma/2)**2/b**2)
        #B =0
        #return epsilon * (A + B)
        return A

    popt, pcov = curve_fit(coeff, x, U, p0 = [8, 8, 5, 2], maxfev=50000)
    print 'sigma= %f p=%f k=%f b=%f'%(popt[0], popt[1],popt[2],popt[3])
    if show:
        xmin=6
        xmax=1.5*popt[0]
        ax1.plot(x,U,'xr')
        ax1.plot(np.linspace(xmin,xmax,50),coeff(np.linspace(xmin,xmax,50),popt[0],popt[1],popt[2],popt[3]),'--r')
    return [popt[0], popt[1],popt[2],popt[3]]

def exp_Potential_fit(ax1,r, U,sigma=20,p=12,epsilon=1,l_c=1,show=False,umin=6):
    def findindex(r):
        for i in range(r.shape[0]):
            if U[i] < umin :
                return i
    def coeff(r, sigma, p, l_c=1):
        espilon=1.0
        B = epsilon*(sigma/r)*np.exp(-(r/sigma)**p)
        return B
    xmin_index= findindex(r)
    xmin_index = 0
    popt, pcov = curve_fit(coeff, r[xmin_index:], U[xmin_index:], p0 = [sigma,
        p], maxfev=100000)
    if show:
        xmax = r[-1]
        xmin = r[0]
        print xmax
        print 'sigma= %f p=%f '%(popt[0], popt[1])
        x = np.linspace(xmin,xmax,100)
        ax1.plot(x,coeff(x,popt[0],popt[1]),'-k',label='fit')
    return [popt[0], popt[1]]

def rexp_Potential_fit(ax1,r, U,sigma=20,p=12,epsilon=1,l_c=1,show=False,umin=6):
    def findindex(r):
        for i in range(r.shape[0]):
            if U[i] < umin :
                return i
    def coeff(r, sigma, p, l_c=1):
        espilon=1.0
        B = epsilon*(sigma/(r-sigma/l_c))*np.exp(-(r/sigma)**p)
        return B
    xmin_index= findindex(r)
    xmin_index = 0
    popt, pcov = curve_fit(coeff, r[xmin_index:], U[xmin_index:], p0 = [sigma,
        p,  l_c], maxfev=100000)
    if show:
        xmax = r[0]
        xmin = r[-1]
        print 'sigma= %f p=%f l_c=%f'%(popt[0], popt[1],popt[2])
        x = np.linspace(xmin,xmax,100)
        ax1.plot(x,coeff(x,popt[0],popt[1],popt[2]),'-k',label='fit')
    return [popt[0], popt[1],popt[2]]

def rexp_Potential_fit_ep(ax1,r, U,sigma,p,epsilon=1,l_c=1,show=False):
    def coeff(r, sigma, p, l_c=1,epsilon=1.0):
        B = epsilon*(sigma/(r-sigma/l_c))*np.exp(-(r/sigma)**p)
        return B
    xmin=6.5
    xmax = 16
    def findindex(r):
        for i in range(r.shape[0]):
            if r[i] > 9.8:
                return i
    xmin_index= findindex(r)
    popt, pcov = curve_fit(coeff, r[xmin_index:], U[xmin_index:], p0 = [sigma,
        p,  l_c, epsilon], maxfev=100000)
    print 'sigma= %f p=%f l_c=%f ep=%f'%(popt[0], popt[1],popt[2],popt[3])
    x = np.linspace(xmin,xmax,100)
    ax1.plot(x,coeff(x,popt[0],popt[1],popt[2],popt[3]),'-g',label='fit_ep')


def Soft_Potential_plotN(ax1,sigma,p,k,b,l_c,show=False):
    epsilon=1
    def coeff(r, sigma, p, k, b):
        #hard core
        c = math.log(1+2**0.5)
        #A =  1 / (np.sinh((c*r/sigma))**p)
        A =  1 / (np.sinh((c*r/8.5))**10)/(c*r/8.5)
        #Soft Core
        # Alex Function
        #w = r / (sigma* l_c)
        #B = k * np.cos(math.pi*((1-b)*w / (b+(1-2*b)*w)))
        #Gumbel Distribution Function
        z = (r - sigma) / l_c
        B = k* np.exp(-(z+np.exp(-z)))
        #B = k*np.exp(-(r-7)**2/b**2) 
    #            (r-sigma)/b**2 *\
    #        np.cos(math.pi*(r-sigma/(sigma*lc)))
        if count == 1:
            return B
        if count == 2:
            return A
        if count == 3:
            return A+B
    #return epsilon * (A + B)
    xmin=0
    xmax=12
    count = 1
    ax1.plot(np.linspace(xmin,xmax,200),coeff(np.linspace(xmin,xmax,200),sigma,p,k,b),
            '--b',label='soft core',lw=2)
    count =2
    ax1.plot(np.linspace(xmin,xmax,200),coeff(np.linspace(xmin,xmax,200),sigma,p,k,b),
            '-.g', label = 'hard core',lw=2)
    count = 3
    ax1.plot(np.linspace(xmin,xmax,200),coeff(np.linspace(xmin,xmax,200),sigma,p,k,b),
            '-k', label = 'H+S',lw=2)

def Soft_Potential_plot_diff(ax1,x,U,sigma,p,k,b):
    epsilon=1
    def coeff(r, sigma, p, k, b):
        c = math.log(1+2**0.5)
        w = 2*r / (sigma* l_c)
        A =  1 / (np.sinh((c*r/sigma))**p)
        #B = k * np.cos(math.pi*((1-b)*w / (b+(1-2*b)*w)))
        z = (r - sigma) / l_c
        B = k* np.exp(-(z+np.exp(-z)))
        return A
    #return epsilon * (A + B)
    v = []
    for i in range(len(x)):
        v.append(U[i] - coeff(x[i],sigma,p,k,b))
    ax1.plot(x,v,'--g')

def Soft_Potential_plot(ax1,sigma,p,epsilon=1,l_c=1,show=False):
    def coeff(r, sigma, p, epsilon=5,l_c=1):
        #A =  epsilon * ( 1 / (np.sinh((r/sigma))**p))/(r/sigma)
        B = epsilon*(sigma/(r-sigma/l_c))*np.exp(-(r/sigma)**p)
        return B
    def coeff_der(r, sigma, p, epsilon=5,l_c=1):
        db_a = -epsilon*sigma*np.exp(-(r/sigma)**p)/(r-sigma/l_c)
        db_b = 1/(r-sigma/l_c)+p*r**(p-1)/sigma**p
        DB = db_a*db_b
        return DB
    xmin=9
    xmax=16
    x = np.linspace(xmin,xmax,50)
    ax1.plot(x,coeff(x,sigma,p,epsilon,l_c),'--b',label='potential-fit')
    #ax1.plot(x,-coeff_der(x,sigma,p,epsilon,l_c),'-b',label='test')
    #U = coeff(x,sigma,p,epsilon,l_c)
    #ax1.plot(x,-np.gradient(U,x[1]-x[0]),'xb',label='test')

if __name__ == '__main__':
    plt.close()
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    fid = open('pmf.inp','r')
    sigma = 10.25
    p=8.32
    k=1
    b=1
    V = []
    F = []
    r = []
    for line in fid.readlines():
        if line[0] != '#':
            x = float(line.split()[0])
            if  x < 14:
                r.append(float(line.split()[0]))
                V.append(float(line.split()[1]))
                F.append(float(line.split()[2]))
    fid.close()
    #Soft_Potential(ax1,np.array(r),np.array(V),show=True)
    #sp=4
    #sigma = 8.68
    #p=8.49
    #epsilon=1.0
    #l_c = 1.281
    sigma = 12
    p=15
    epsilon=.6
    l_c = 1.32
    Soft_Potential_plot(ax1,sigma,p,epsilon,l_c)
    rexp_Potential_fit_ep(ax1,np.array(r),np.array(V),sigma,p,epsilon,l_c)
    rexp_Potential_fit(ax1,np.array(r),np.array(V),sigma,p,epsilon,l_c)
    ax1.plot(r,V,'xr',label="pmf")
    #ax1.plot(r,-np.gradient(V,r[1]-r[0]),label='gradient')
    #ax1.plot(r,F,'ob',label="Force")
    ax1.set_yscale('log')
    #ax1.set_xscale('log')
    ax1.set_ylim((.01,100))
    ax1.set_xlim((9,14))
    plt.legend(loc=1)
    ax1.set_xlabel('r')
    ax1.set_ylabel('U(r)')
    plt.show()
