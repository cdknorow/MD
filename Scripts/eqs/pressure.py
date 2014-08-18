import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MD

def pressure():
    file1 = 'pressure.log'
    file2 = 'box_length.txt'
    fid = open(file1,'r')
    pressure = [[] for i in range(7)]
    fid.readline()
    Ptotal = []
    for P in fid.readlines():
        for i in range(7):
            pressure[i].append(float(P.split()[i+1]))
        Ptotal.append((pressure[1][-1]+pressure[2][-1]+pressure[3][-1])/3.)
    print pressure
    fid.close()
    fid = open(file2,'r')
    box_length = []
    fid.readline()
    for L in fid.readlines():
        box_length.append(float(L.split()[1]))
    fid.close()

    #plot the pressures vs time
    x = range(len(pressure[1]))
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, pressure[0], label='pressure')
    ax1.plot(x, pressure[1], label='xx')
    ax1.plot(x, pressure[2], label='yy')
    ax1.plot(x, pressure[3], label='zz')
    ax1.plot(x, Ptotal, label='psum/3')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Pressure')
    ax1.legend(loc=4)
    plt.savefig('pressures.png')
    #plt.show()
    plt.close()
    #plot the pressure vs box size
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(box_length,pressure[0],'x',label='Pressure')
    ax2.set_xlabel('Volume')
    ax2.set_ylabel('Pressure')
    ax2.legend(loc=4)
    plt.savefig('PvsV.png')
    plt.show()
    plt.close()
    x = range(len(pressure[1]))
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x, pressure[0], label='pressure')
    ax1.plot(x, pressure[4], label='xy')
    ax1.plot(x, pressure[5], label='yz')
    ax1.plot(x, pressure[6], label='xz')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Pressure')
    ax1.legend(loc=4)
    plt.savefig('pressuresxy.png')
    plt.show()
    plt.close()

def print_box_volume(L_cont,delta=10):
    fid = open('box_length.txt','w')
    fid.write('frame  L\n') 
    for i in range(0,len(L_cont)):
        fid.write(('%i %.2f %.2f %.2f\n')%(i,L_cont[i][0],L_cont[i][1],L_cont[i][2]))
    fid.close()

if __name__ == "__main__":
    M=MD.ReadCord()
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    print_box_volume(L_cont)
    pressure()
