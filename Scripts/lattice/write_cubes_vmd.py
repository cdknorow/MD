# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import math
sys.path.append('/home/cdknorow/Dropbox/Software/')
import MD
import MD.util as util
import MD.base.points as points
import draw_cube as dw
#allign the vectors so that they are drawn in the correct direction
def allign(v1,v2,v3):
    vc = [v1,v2,v3]
    us = np.array([1,1,1])
    vs = np.array([1,1,-1])
    ws = np.array([1,-1,1])
    usmall = 0
    vsmall = 0
    wsmall = 0
    angle = []
    for j in range(3):
        angle.append([points.angle_between(us,vc[j]),j])
    angle.sort()
    if angle[0][0] <1.0:
        print 'swapping'
        if points.angle_between(us,-vc[angle[-1][1]])<1.0:
            vc[angle[-1][1]] = -vc[angle[-1][1]]
            print points.angle_between(us,vc[angle[-1][1]])
        if points.angle_between(us,-vc[angle[-2][1]])<1.0:
            vc[angle[-2][1]] = -vc[angle[-2][1]]
            print points.angle_between(us,vc[angle[-2][1]])
    x = [0,1,2]
    for j in x:
        if points.angle_between(us,vc[j])<=points.angle_between(us,vc[usmall]):
            usmall = j
    del x[usmall]
    wsmall = x[0]
    print x
    for j in x:
        if points.angle_between(ws,vc[j])<=points.angle_between(ws,vc[wsmall]):
            wsmall = j
    del x[x.index(wsmall)]
    vsmall = x[0]
    print usmall,wsmall,vsmall
    return vc[usmall], vc[wsmall], vc[vsmall]
def read_xyz():
    fid = open('dna.xyz','r')
    N_atoms = fid.readline()
    fid.readline()
    M = fid.readlines()
    NP_center = []
    NP_name = []
    Z = []
    #NP_center = index of center of particle
    #V = numpy array of all particle positions
    index = []
    for line in M:
        p = np.array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
        if line.split()[0] == 'V' or line.split()[0] == 'W':
            NP_center.append(p)
            NP_name.append(line.split()[0])
        if line.split()[0] == 'Z':
            Z.append(p)
    return NP_center,Z,NP_name
def read_color_index(frame = ''):
    if frame == '':
        fid = open('animate/bakos_index.txt','r')
    else:
        fid = open('animate/bakos_index%i.txt'%frame,'r')
    M = fid.readlines()
    A = []
    for line in M:
        A.append(float(line))
    return A
def run_colored(frame=1):
    ##FROM DCD
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    #last = M.frames
    #V_index = M.get_index(['V'])
    #W_index = M.get_index(['W'])
    ## 
    #Z = util.pickle_load('Z.pkl')
    #VW = util.pickle_load('VW.pkl')
    ##
    #k = last
    #FROM XYZ
    VW,Z,VW_name = read_xyz()
    fid = open('square%i.tcl'%frame,'w')
    #set x to be the particles you wish to choose
    #x = np.array([47190, 77220, 72215, 80795, 42185, 17875, 72930, 52910,
    #    67925])
    #index 47190 77220 72215 80795 42185 17875 72930 52910 67925
    #x = x/715
    #print x
    #x = [1,2,3,4,5]
    #for i in x:
    A = read_color_index()
    for i in range(len(VW)):
        c = VW[i]
        v1 = Z[i*6] - c
        v2 = Z[i*6+2] - c
        v3 = Z[i*6+4] - c
        #v1 = Z[i*6] - c
        #v2 = Z[i*6+1] - c
        #v3 = Z[i*6+5] - c
        #v1,v2,v3 = allign(v1,v2,v3)
        color = ['blue','red','green','orange','purple']
        dw.draw_cube(fid,c,v1,v2,v3,color[A[i]])
    fid.close()
def animate():
    ##FROM DCD
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    Z = util.pickle_load('Zbound.pkl')
    VW = util.pickle_load('VW.pkl')
    animate_script = open('animate/animate.tcl','w')
    for i in range(VW.shape[1]):
        if i < VW.shape[1]/2:
            A.append(0)
        else:
            A.append(1)

    #for k in range(len(VW.shape[0])):
    for k in range(50,100):
        fid = open('animate/square%i.tcl'%k,'w')
        #A = read_color_index(frame=k)
        for i in range(VW.shape[1]):
            c = VW[k][i]
            v1 = Z[k][i*6] - c
            v2 = Z[k][i*6+2] - c
            v3 = Z[k][i*6+4] - c
            #v1 = Z[k][i*6] - c
            #v2 = Z[k][i*6+1] - c
            #v3 = Z[k][i*6+5] - c
            #v1,v2,v3 = allign(v1,v2,v3)
            color = ['blue','red','green','orange','purple']
            dw.draw_cube(fid,c,v1,v2,v3,color[int(A[i])])
        animate_script.write('\nsource square%i.tcl\n display update ui\n'%k)
        animate_script.write('\nafter 1000\ndraw delete all\n display update ui\n')
        fid.close()

#For single directories
if __name__ == '__main__':
    #run_colored()
    animate()

