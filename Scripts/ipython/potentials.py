import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
matplotlib.rc('font', **{'size':'20'})
plt.close()
xlabel = 'r'
ylabel = 'U(r)'
sigma=1.0
epsilon=1.0
xmin = 0.7
Uwca = []
rwca = np.linspace(xmin,2**(1./6),num=100)
for r in rwca:
    Uwca.append( 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6 + 1/4.))
    
U20 = []
U12 = []
U8 = []
U6 = []
U6M = []
rpower = np.linspace(xmin,5,num=3000)
for r in rpower:
	U20.append( 4 * epsilon * ( (sigma / r)**20))
	U12.append( 4 * epsilon * ( (sigma / r)**12))
	U8.append( 4 * epsilon * ( (sigma / r)**8))
	U6.append( 4 * epsilon * ( (sigma / r)**6))
	U6M.append( epsilon * ( (sigma / r)**6))

#for r in rpower:
#    U20.append( math.exp(-4 * epsilon * ( (sigma / r)**20)))
#    U12.append( math.exp(-4 * epsilon * ( (sigma / r)**12)))
#    U8.append( math.exp(-4 * epsilon * ( (sigma / r)**8)))
#    U6.append( math.exp(-4 * epsilon * ( (sigma / r)**6)))
#    U6M.append( math.exp(-epsilon * ( (sigma / r)**6)))
sigma = 2.17
p = 12.3604
kappa = .1
lp = 0.7
k=10

def Usoft(r, p, sigma, epsilon=1.0):
    V = epsilon * ( 1 / (np.sinh((r/sigma))**p))
    return (V)
 
def Uharm(r, sigma,lp,k,epsilon=1.0):
    V = epsilon * k * ((r/sigma - lp)/lp)**2
    return (V)

print sigma*lp
rsoft = np.linspace(1,5,num=300)
rharm = np.linspace(1,2*sigma*lp,num=300)
fig1 = plt.figure(figsize=(15,7), dpi=100)
ax1 = fig1.add_subplot(111)
ax1.plot(rsoft,Usoft(rsoft,p=p,sigma=sigma),label='Usoft')
for i in range(5):
	print k*i
	ax1.plot(rharm,Uharm(rharm,sigma=sigma*2,lp=lp,k=k*i),label='Uharm%i'%k*i)
ax1.set_ylim(0,5)
ax1.set_xlim(1.5,5)
ax1.legend(loc=1)
plt.show()

if False:
	fig1 = plt.figure(figsize=(15,7), dpi=100)
	ax1 = fig1.add_subplot(111)
	ax1.plot(rpower,U20,label='U20')
	ax1.plot(rpower,U12,label='U12')
	ax1.plot(rpower,U8,label='U8')
	ax1.plot(rpower,U6,label='U6')
	ax1.plot(rpower,U6M,label='U6M')
	ax1.plot(rwca,Uwca,label='wca')
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(ylabel)
	ax1.set_ylim(0,.02)
	ax1.set_xlim(.7,5)

	ax1.legend(loc=1)
	plt.show()