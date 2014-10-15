
import numpy as np
import os
import MD.base.points as points



## \brief get the average position of particles taking periodic boundaries into account
#
# \return average position of particles over time
#
# \param M a group of atoms in matrix form [frames][atoms][x,y,z]

def average_position(VW, n_start, n_frames, L, GenLattice = False):
    if n_frames == 1:
        return VW[n_start]
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

    if GenLattice:
        fid = open('average.xyz','w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
           fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()

    return Average


## \brief get the average position of particles don't accouunt for boundary conditions
#
# \return average position of particles over time
#
# \param M a group of atoms in matrix form [frames][atoms][x,y,z]

def average_position_simple(M):
    A = np.zeros((M.shape[1],3))
    for i in range(M.shape[1]):
        count = 1
        p = np.array([0.0,0.0,0.0])
        p += M[0][i]
        for j in range(1,M.shape[0]):
            if points.dist_np(M[0][i],M[j][i])[0] <  5:
                p += M[j][i]
                count += 1
        p = p/count
        A[i] = p
    return A
