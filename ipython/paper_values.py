import numpy as np
from scipy import integrate as integrate 
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font', **{'size':'20'})
fig1 = plt.figure(figsize=(12,8), dpi=100)
ax1 = fig1.add_subplot(111)
gamma_potential = np.array([0.1,0.25,0.5,1.0,1.5])
#N=256
#(PE-PE_0)/NkT
potential_paper = np.array([0.247,0.577,1.036,1.739,2.195])
Z_paper = []
for i in range(potential_paper.shape[0]):
	Z_paper.append(2*(potential_paper[i] + 7.227*gamma_potential[i]**2))
F_potential_paper = integrate.cumtrapz(Z_paper, gamma_potential,initial=0)
ax1.plot(gamma_potential,F_potential_paper, 'o',label='Calculated Free Energy')

real_data_gamma = [0.3162,0.6324,0.9487,1.2649,1.5811,1.8974,2.2136,2.317]
real_data_FE = [1.4609,3.5364,6.2504,9.6153,13.642,18.3382,23.706,26.644]
ax1.plot(real_data_gamma,real_data_FE, 'o',label='Paper Free Energy')
plt.show()