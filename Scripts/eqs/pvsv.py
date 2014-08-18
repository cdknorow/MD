import sys
import os
import numpy as np
from scipy import integrate as integrate 
import matplotlib.pyplot as plt


Dirlist = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
Dirlist.sort()
Dlist = Dirlist

file1 = 'pressure.log'
file2 = 'box_length.txt'
box_length = [[] for i in range(len(Dlist))]
pressure = [[[] for i in range(7)] for i in range(len(Dlist))]
for i in range(len(Dlist)):
    os.chdir(Dlist[i])
    fid = open(file1,'r')
    fid.readline()
    fid.readline()
    count = 0
    for P in fid.readlines():
        for j in range(7):
            if count ==0:
                pressure[i][j].append(float(P.split()[1+j]))
        count+=1
        if count==5:
            count=0
    fid.close()
    fid = open(file2,'r')
    fid.readline()
    fid.readline()
    for L in fid.readlines():
        box_length[i].append(float(L.split()[1]))
    fid.close()
    if os.path.isdir('cont'):
        os.chdir('cont')
        print os.getcwd()
        fid = open(file1,'r')
        fid.readline()
        fid.readline()
        for P in fid.readlines():
            for j in range(7):
                pressure[i][j].append(float(P.split()[1+j]))
        fid.close()
        fid = open(file2,'r')
        fid.readline()
        fid.readline()
        for L in fid.readlines():
            box_length[i].append(float(L.split()[1]))
        fid.close()
        if os.path.isdir('cont'):
            os.chdir('cont')
            print os.getcwd()
            fid = open(file1,'r')
            fid.readline()
            fid.readline()
            for P in fid.readlines():
                for j in range(7):
                    pressure[i][j].append(float(P.split()[1+j]))
            fid.close()
            fid = open(file2,'r')
            fid.readline()
            fid.readline()
            for L in fid.readlines():
                box_length[i].append(float(L.split()[1]))
            fid.close()
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')


Pavg = [[] for i in range(len(Dlist))]
Lavg = [[] for i in range(len(Dlist))]
Perror = [[] for i in range(len(Dlist))]
for i in range(len(Dlist)):
    if len(pressure[i][0])!=len(box_length[i]):
        del pressure[i][0][-1]
        print len(pressure[i][0])
        print len(box_length[i])
for i in range(len(Dlist)):
    L = box_length[i][0]
    psum = 0
    count = 0
    for j in range(len(box_length[i])-1):
        if box_length[i][j] == L:
            psum += pressure[i][0][j]
            count+=1
        else:
            if L <= 92:
                Pavg[i].append(psum/count)
                Lavg[i].append(L)
                #calculate the error
                counter = 0
                perr = 0
                for k in range(len(pressure[i][0])):
                    try:
                        if box_length[i][k] == L:
                            perr += (Pavg[i][-1]-pressure[i][0][k])**2
                            counter += 1
                    except:
                        print 'box length and pressure length are not the same'
                        print 'box length'
                        print len(box_length[i])
                        print 'pressure length'
                        print len(pressure[i][0])
                Perror[i].append((perr/counter)**0.5)
            psum = pressure[i][0][j]
            L = box_length[i][j]
            count = 1
#sp8
print Dlist
if Dlist[0].split('_')[0] == 'sp8':
    Dlabel =   ['tric_1e6','bcc_1e6','tric_2e6','bcc_5e6',
                 'tric_5e6','tric_5e6','tric_5e6','bcc_5e6']
    Dmarker = ['ro','gs','b^','b+','k,','y1','r2','g3','o4']
#sp4
else:
    Dlabel =['bcc_5e6','tric_1e6','tric_2.5e6','tric_1e6','tric_5e6',
              'bcc_5e6','tric_5e6','tric_1e6']
    Dmarker = ['ro','gs','b^','b+','k,','y1','r2','g3','o4']
#plot the pressure vs box size

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for i in range(len(Dlist)):
    #ax2.errorbar(Lavg[i],Pavg[i],yerr=Perror[i], fmt=Dmarker[i],label=Dlabel[i])
    #if i ==0:
    ax2.scatter(box_length[i],pressure[i][0],label=Dlabel[i])
    #else:
    #    ax2.scatter(box_length[i],pressure[i][0],label=Dlabel[i])
ax2.set_xlabel('Volume')
ax2.set_ylabel('Pressure')
ax2.legend(loc=1)
plt.savefig('PvsVavg.png')
#plt.ylim(0,0.7)
#plt.xlim(87,91)
plt.show()
plt.close()
#plot the change in free energy
#DF = -Int^Vf_Vi{PdV}
Favg = [[] for i in range(len(Dlist))]
Vavg = [[] for i in range(len(Dlist))]
Vmid = [[] for i in range(len(Dlist))]
for i in range(len(Pavg)):
    for j in range(len(Lavg[i])):
        Vavg[i].append(Lavg[i][j]**3)
        if j != 0:
            Vmid[i].append((Vavg[i][j]+Vavg[i][j-1])/2)
    Favg[i]= -integrate.cumtrapz(Pavg[i],Vavg[i])
    print Favg[i]
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for i in range(len(Dlist)):
    ax2.plot(Vmid[i],Favg[i],'-x',label=Dlabel[i])
ax2.set_xlabel('Volume')
ax2.set_ylabel('Free Energy')
ax2.legend(loc=1)
plt.savefig('FvsVavg.png')
#plt.ylim(0,0.7)
#plt.xlim(87,91)
plt.show()
plt.close()
