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
############ 
# This will shift particle positions in a way that removes the periodic boundary
# effects
#########
reload(MD)
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
##\brief unwrap particles with bonds from periodic boundaries, creating a system
# entirley contained within the box 
def remove_periodic(frame=1):
    #print out directory
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = M.box_volume()
    print L[0]
    last = M.frames
    print M.num_particles
    delta = 1
    V=M.all_cord_range(start=1,delta=1,last=2)
    print V
    NP_center = []
    for i,j in enumerate(M.universe):
        if M.universe[i] == "V" or M.universe[i] == 'W':
            NP_center.append(i)
    R = copy.deepcopy(V)
    #the center of the partilcl
    C = NP_center[0]
    NP_size = NP_center[1] - NP_center[0]
    NP2_size = NP_center[-1] - NP_center[-2]
    for count,k in enumerate(NP_center):
        #the number of particles in a NP
        if count<(len(NP_center)/2):
            print k
            for i in range(NP_size):
                #apply boundary conditions to move particle
                print k
                V[0][i+k-C] = boundary(V[0][i+k-C],V[0][k],L[0])
        else:
            for i in range(NP2_size):
                #apply boundary conditions to move particle
                V[0][i+k-C] = boundary(V[0][i+k-C],V[0][k],L[0])



    print 'writing dna.xyz'
    out = open('dna%i.xyz'%frame,'w')
    out.write(('%i\n\n')%(V.shape[1]))
    count = 0
    for i in range(V.shape[1]):
        out.write(('%s %.2f %.2f %.2f\n')%(M.universe[count],V[0][count][0],V[0][count][1],V[0][count][2]))
        count+=1

    #for i in range(V.shape[1]):
    #    x = V[0][i][1]
    #    y = -V[0][i][0]
    #    z = V[0][i][2]

    #    V[0][i][0] = x
    #    V[0][i][1] = y
    #    V[0][i][2] = z

    #print 'writing dna.xyz'
    #out = open('dna_r.xyz','w')
    #out.write(('%i\n\n')%(V.shape[1]))
    #count = 0
    #for i in range(V.shape[1]/2):
    #    if M.universe[count] == 'N':
    #        out.write(('Q %.2f %.2f %.2f\n')%(V[0][count][0],V[0][count][1],V[0][count][2]))
    #    else:
    #        out.write(('%s %.2f %.2f %.2f\n')%(M.universe[count],V[0][count][0],V[0][count][1],V[0][count][2]))
    #    count+=1
    #for i in range(V.shape[1]/2,V.shape[1]):
    #    if M.universe[count] == 'N':
    #        out.write(('N %.2f %.2f %.2f\n')%(V[0][count][0],V[0][count][1],V[0][count][2]))
    #    else:
    #        out.write(('%s %.2f %.2f %.2f\n')%(M.universe[count],V[0][count][0],V[0][count][1],V[0][count][2]))
    #    count+=1
##\brief remove a particle that is too close to a neighboring particle 
#after things are merged via packmol
def nearpart_remove():
    #print out directory
    fid = open('dna_2nuc.xyz','r')
    N_atoms = fid.readline()
    fid.readline()
    M = fid.readlines()
    V = []
    count = 0
    NP_center = []
    #NP_center = index of center of particle
    #V = numpy array of all particle positions
    #index = "name of particle at each index'
    index = []
    for line in M:
        index.append(line.split()[0])
        V.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
        if line.split()[0] == 'V' or line.split()[0] == 'W':
            NP_center.append(count)
        count += 1
    fid.close()
    V = np.array(V)
    Lx = 220
    Ly = 110
    Lz = 110
    L = np.array([Lx, Ly, Lz])
    print L
    #the center of the partilcl
    too_close = []
    print NP_center
    C = NP_center[0]
    NP_size = NP_center[1] - NP_center[0]
    for i in NP_center:
        print i/NP_size
        for j in NP_center:
            if i != j:
                if points.dist(V[i],V[j],L)[0]<12:
                    too_close.append([min(i,j),max(i,j)])
    too_close = points.unique(too_close)
    remove = []
    for i in too_close:
        remove.append(i[0])

    #remove the NP which are too_close!
    #to do and write the new dna.xyz file
    print 'writing dna.xyz'
    out = open('dna_close_remove.xyz','w')
    out.write(('%i\n\n')%(V.shape[0]))
    count = 0
    for i in range(V.shape[0]):
        if i in remove:
            i+=1
        if count>=V.shape[0]-1:
            print 'breaking'
            break
        for j in range(NP_size):
            count = i * NP_size + j
            out.write(('%s %.2f %.2f %.2f\n')%(index[count],V[count][0],V[count][1],V[count][2]))
##\brief shift particles so they are in a new box, causing all particles outside of the
#of certain distance to wrap around the box
def shift_boundary():
    #print out directory
    fid = open('dna_close_remove.xyz','r')
    N_atoms = fid.readline()
    fid.readline()
    M = fid.readlines()
    V = []
    count = 0
    NP_center = []
    #NP_center = index of center of particle
    #V = numpy array of all particle positions
    #index = "name of particle at each index'
    index = []
    out_x = 0
    out_y = 0
    out_z = 0
    for line in M:
        index.append(line.split()[0])
        V.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
        if line.split()[0] == 'V' or line.split()[0] == 'W':
            NP_center.append(count)
        if line.split()[0] == 'N':
            print V[-1]
            if out_x < abs(V[-1][0]):
                out_x = abs(V[-1][0])
            if out_y < abs(V[-1][1]):
                out_y = abs(V[-1][1])
            if out_z < abs(V[-1][2]):
                out_z = abs(V[-1][2])
        count += 1
    fid.close()
    V = np.array(V)
    Lx = out_x*2+8
    Ly = out_y*2+8
    Lz = out_z*2+8
    L = np.array([Lx, Ly, Lz])
    print L
    #todo, force particles outside of distance L to wrap around box
    for i in range(V.shape[0]):
        #apply boundary conditions to move particle
        V[i] = flip_boundary(V[i],L)


    print 'writing dna.xyz'
    out = open('dna_boundary_shrink.xyz','w')
    out.write(('%i\n\n')%(V.shape[0]))
    count = 0
    for i in range(V.shape[0]):
        out.write(('%s %.2f %.2f %.2f\n')%(index[i],V[i][0],V[i][1],V[i][2]))

if __name__ == '__main__':
    #run_shift()
    #run_nearpart_remove()
    #shift_boundary()
    remove_periodic()

