
## \pacakge MD.dna_scripts.position_shift
#
# \brief this script manipulates the positions of particles with respect to periodic
# boundaries
import sys
sys.path.append('/home/cdknorow/Dropbox/Software')
import math
import os
import numpy as np
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import copy
import draw_cube as dw
import MD.analysis.cubatic_order as cub
############ 
# This will shift particle positions in a way that removes the periodic boundary
# effects
#########
reload(MD)
# M [cube][x,y,z]
# Z [edge][x,y,z], edge 1-6 represent first cube
def Average(VW,L,n_start,n_frames,write=True,save='star.xyz'):
    #Average_position
    Average = np.zeros(VW[0].shape)
    for k in range(n_start,n_start+n_frames):
        for i in range(VW.shape[1]):
            #check the distance between
            for j in range(3):
                if abs(VW[k][i][j]-VW[n_start][i][j]) < L[k][j]/2:
                    Average[i][j] += VW[k][i][j]
                elif VW[k][i][j]-VW[n_start][i][j] < -L[k][j]/2:
                    Average[i][j] += VW[k][i][j]+L[k][j]
                elif VW[k][i][j]-VW[n_start][i][j] > L[k][j]/2:
                    Average[i][j] += VW[k][i][j]-L[k][j]
   # fix any points that may be outside of the box after averaging
    for i in range(Average.shape[0]):
        for j in range(3):
            Average[i][j] /= n_frames
            if Average[i][j] > L[n_start][j]:
                Average[i][j] -= L[n_start][j]
            if Average[i][j] < -L[n_start][j]:
                Average[i][j] += L[n_start][j]
    if write:
        fid = open(save,'w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()
    return Average
def cube_vectors(V,E,L):
    vec = []
    for i in E:
        vec.append(points.unit(points.vector1d(V,i,L)))
    return vec
def cubic_order_animate(VW,Z,L_cont,x):
    Q = []
    reload(cub)
    xvalue = []
    bakos_pick=[]
    frame = 40
    bakos_pick = [161,476]
    bakos_cubes = []
    for index in bakos_pick:
        bakos_cubes.append(cube_vectors(VW[frame][index],Z[frame][index*6:index*6+6],L_cont[frame]))
    for i,k in enumerate(x):
        print i,k
        cub.least_fit_cubes(VW[k],Z[k],L_cont[k],k,bakos_pick=bakos_cubes)
        #cub.least_fit_cubes(VW[k],Z[k],L_cont[k],k,bakos_pick=True)
#check boundary conditions
def boundary(a,b,L):
    if a[0] - b[0] > L[0]/2:
        a[0] = a[0] - L[0]

    if a[0] - b[0] < -L[0]/2:
        a[0] = a[0] + L[0]

    if a[1] - b[1] > L[1]/2:
        a[1] = a[1] - L[1]

    if a[1] - b[1] < -L[1]/2:
        a[1] = a[1] + L[1]

    if a[2] - b[2] > L[2]/2:
        a[2] = a[2] - L[2]

    if a[2] - b[2] < -L[2]/2:
        a[2] = a[2] + L[2]

    return a
#check boundary conditions
def flip_boundary(a,L):
    if a[0] > L[0]/2:
        a[0] = a[0] - L[0]

    if a[0] < -L[0]/2:
        a[0] = a[0] + L[0]

    if a[1] > L[1]/2:
        a[1] = a[1] - L[1]

    if a[1] < -L[1]/2:
        a[1] = a[1] + L[1]

    if a[2] > L[2]/2:
        a[2] = a[2] - L[2]

    if a[2] < -L[2]/2:
        a[2] = a[2] + L[2]

    return a
# -*- coding: utf-8 -*-
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
def animate(frames=[0,1]):
    #print out directory
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = M.box_volume()
    last = M.frames
    try:
        Z = util.pickle_load('Z.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
        Z=M.cord_auto(['Z'])
        VW=M.cord_auto(['V','W'])
        util.pickle_dump(Z,'Z.pkl')
        util.pickle_dump(VW,'VW.pkl')
    cubic_order_animate(VW,Z,L,frames)
    for k in frames:
    #apply boundary conditions to move particle
        for i in range(VW.shape[1]):
            for j in range(i*6,i*6+6):
                Z[k][j] = boundary(Z[k][j],VW[k][i],L[k])
    util.pickle_dump(Z,'Zbound.pkl')
    animate_script = open('animate/animate.tcl','w')
    #for binary systems use this A
    #A = []
    #for i in range(VW.shape[1]):
    #    if i < VW.shape[1]/2:
    #        A.append(0)
    #    else:
    #        A.append(1)
    #for k in range(len(VW.shape[0])):
    for k in frames:
        fid = open('animate/square%i.tcl'%k,'w')
        fid2= open('animate/square%i.txt'%k,'w')
        A = read_color_index(frame=k)
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
            dw.draw_cube(fid,fid2,c,v1,v2,v3,color[int(A[i])])
        cub.gaussmap_color(VW,Z,L,k,A)
        animate_script.write('\nsource square%i.tcl\n display update ui\n'%k)
        animate_script.write('\nafter 1000\ndraw delete all\n display update ui\n')
        fid.close()
def animate_binary(frames=[0,1]):
    #print out directory
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = M.box_volume()
    last = M.frames
    try:
        Z = util.pickle_load('Z.pkl')
        VW = util.pickle_load('VW.pkl')
        npType = M.get_type(['V','W'])
    except:
        Z=M.cord_auto(['Z'])
        VW=M.cord_auto(['V','W'])
        npType = M.get_type(['V','W'])
        util.pickle_dump(Z,'Z.pkl')
        util.pickle_dump(VW,'VW.pkl')
    for k in frames:
    #apply boundary conditions to move particle
        for i in range(VW.shape[1]):
            for j in range(i*6,i*6+6):
                Z[k][j] = boundary(Z[k][j],VW[k][i],L[k])
    util.pickle_dump(Z,'Zbound.pkl')
    animate_script = open('animate/animate.tcl','w')
    #for binary systems use this A
    #fid = open('sc_num.txt','r')
    #crystal = []
    #for line in fid.readlines():
    #    crystal.append(float(line.split()[0]))
    #fid.close()
    A = []
    for i in range(VW.shape[1]):
        if npType[i] == 'V':
            A.append(0)
        else:
            A.append(1)
    for k in frames:
        fid = open('animate/square%i.tcl'%k,'w')
        for i in range(VW.shape[1]):
            c = VW[k][i]
            #v1 = Z[k][i*6] - c
            #v2 = Z[k][i*6+2] - c
            #v3 = Z[k][i*6+4] - c
            v1 = Z[k][i*6] - c
            v2 = Z[k][i*6+1] - c
            v3 = Z[k][i*6+5] - c
            #v1,v2,v3 = allign(v1,v2,v3)
            color = ['blue','red','green','orange','purple']
            dw.draw_cube(fid,c,v1,v2,v3,color[int(A[i])])
        cub.gaussmap_color(VW,Z,L,k,A)
        animate_script.write('\nsource square%i.tcl\n display update ui\n'%k)
        animate_script.write('\nafter 1000\ndraw delete all\n display update ui\n')
        fid.close()
def animate_binary_D(frames=[0,1]):
    #print out directory
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = M.box_volume()
    last = M.frames
    try:
        Z = util.pickle_load('Z.pkl')
        VW = util.pickle_load('VW.pkl')
        npType = M.get_type(['V','W'])
    except:
        Z=M.cord_auto(['Z'])
        VW=M.cord_auto(['V','W'])
        npType = M.get_type(['V','W'])
        util.pickle_dump(Z,'Z.pkl')
        util.pickle_dump(VW,'VW.pkl')
    for k in frames:
    #apply boundary conditions to move particle
        for i in range(VW.shape[1]):
            for j in range(i*6,i*6+6):
                Z[k][j] = boundary(Z[k][j],VW[k][i],L[k])
    util.pickle_dump(Z,'Zbound.pkl')
    animate_script = open('animate/animate.tcl','w')
    #for binary systems use this A
    fid = open('sc_num.txt','r')
    crystal = []
    for line in fid.readlines():
        crystal.append(float(line.split()[0]))
    fid.close()
    A = []
    for i in range(VW.shape[1]):
        if npType[i] == 'V':
            A.append(0)
        else:
            A.append(1)
    for k in frames:
        fid = open('animate/Asquare%i.tcl'%k,'w')
        for i in range(VW.shape[1]):
            c = VW[k][i]
            #v1 = Z[k][i*6] - c
            #v2 = Z[k][i*6+2] - c
            #v3 = Z[k][i*6+4] - c
            v1 = Z[k][i*6] - c
            v2 = Z[k][i*6+1] - c
            v3 = Z[k][i*6+5] - c
            #v1,v2,v3 = allign(v1,v2,v3)
            color = ['blue','red','green','orange','purple']
            if i in crystal:
                dw.draw_cube(fid,c,v1,v2,v3,color[int(A[i])],material='AOShiny')
        cub.gaussmap_color(VW,Z,L,k,A)
        animate_script.write('\nsource square%i.tcl\n display update ui\n'%k)
        animate_script.write('\nafter 1000\ndraw delete all\n display update ui\n')
        fid.close()
    for k in frames:
        fid = open('animate/Bsquare%i.tcl'%k,'w')
        for i in range(VW.shape[1]):
            c = VW[k][i]
            #v1 = Z[k][i*6] - c
            #v2 = Z[k][i*6+2] - c
            #v3 = Z[k][i*6+4] - c
            v1 = Z[k][i*6] - c
            v2 = Z[k][i*6+1] - c
            v3 = Z[k][i*6+5] - c
            #v1,v2,v3 = allign(v1,v2,v3)
            color = ['blue','red','green','orange','purple']
            if i not in crystal:
                dw.draw_cube(fid,c,v1,v2,v3,color[int(A[i])],material='Edgy')
        cub.gaussmap_color(VW,Z,L,k,A)
        animate_script.write('\nsource square%i.tcl\n display update ui\n'%k)
        animate_script.write('\nafter 1000\ndraw delete all\n display update ui\n')
        fid.close()
def animate_beads(frames=[0,1]):
    #print out directory
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = M.box_volume()
    last = M.frames
    try:
        VW = util.pickle_load('VW.pkl')
        npType = M.get_type(['V','W'])
    except:
        VW=M.cord_auto(['V','W'])
        npType = M.get_type(['V','W'])
        util.pickle_dump(VW,'VW.pkl')
    #for binary systems use this A
    for k in frames:
        fid = open('animate/Color%i.tcl'%k,'w')
        fid.write('%i\n\n'%VW.shape[1])
        V = Average(VW,L,k-5,5)
        A = read_color_index(frame=k)
        for i in range(VW.shape[1]):
            x = V[i][0]
            y = V[i][1]
            z = V[i][2]
            color = ['B','R','G','O','P']
            fid.write('%c %.2f %.2f %.2f\n'%(color[int(A[i])],x,y,z))
        fid.close()
if __name__ == '__main__':
    #run_shift()
    #run_nearpart_remove()
    #shift_boundary()
    frames = [40]
    print 'remove priodic'
    print 'create inex map'
    print 'animating'
    animate(frames)
    #animate_beads(frames)
    #animate_binary(frames)
    print 'all done!'

