
## \package MD.analysis.msd 
# \brief This module msd of atoms and takes into account drift of system by
# subtracting out center of mass.
#

import numpy as np
import MD.base.points as points
from MD.util.util import acypher


## \brief Find the distance and direction of drift for a set of particles away
#         from the center point.
#
# This function is used to find the drift of particles from the center of 
# mass position throughout a simulation run. This information is typically
# used by msd to subtract out any artificial movement that occurs during
# a simulation run
# 
# \returns Drift in x,y,z and Drift total
#  
# \param A matrix of [frames][atoms][x,y,z]
# \param L length of box
#
def drift_position(A,L):
    num = A[0].shape[0]
    print "Calculating drift of position"
    D = np.zeros(((A.shape[0]-1),3))
    Dsum = np.zeros(((A.shape[0]-1),1))
    x_whole = 0; y_whole = 0; z_whole = 0
    for k in range(A.shape[0] - 1):
        x_sum = 0; y_sum = 0; z_sum = 0
        for i in range(num):
            d, dx = points.dist(A[k][i], A[k+1][i], L)
            x_sum += dx[0]
            y_sum += dx[1]
            z_sum += dx[2]
        D[k][0] = (x_sum / num + x_whole)
        D[k][1] = (y_sum / num + y_whole)
        D[k][2] = (z_sum / num + z_whole)
        Dsum[k] = (((x_sum / num +x_whole)**2
                    +(y_sum / num + y_whole)**2
                    +(z_sum / num + z_whole)**2)**0.5)
        x_whole = x_sum / num + x_whole
        y_whole = y_sum / num + y_whole
        z_whole = z_sum / num + z_whole
    return  D, Dsum


## \internal
# \brief Find if the particle has crossed the boundary and move it back to its
# place
#
# \returns updated information about crossing the boundary
#
# \param p1 point one
# \param p2 point two
# \param B stores information for how many times boundary has been crossed
#
def boundary(x,L):
    #The particle could be very far past the boundary so 
    #we need to make sure it is put in the right place
    while x > L / 2.0:
        x = x - L
    while x <= -L / 2.0:
        x = x + L
    return x
## \brief Find the msd for a set of particles away from the center point
#
# msd is defined as 
# \f{eqnarray*}
# MSD(t) = \sum_{i=0}^N(r(0)_{i}-r(t)_i)^{0.5} 
# \f}
#where r0 is the inital
# r vector of the particle and ri is ith particle in frame k. the msd is
# computed for each frame .
#
#the drift of the system is also subtracted away to get an accurate measure 
# of the movement of the particles
# 
# \returns time steps and corresponding MSD values
#
# \param A matrix of [frames][atoms][x,y,z]  
# \param L length of box
# \param step time step between dump files
#
def eliminate_drift(A,L,index=False,step=1):
    #find the drift
    Drift = drift_position(A,L)[0]
    #find msd
    num = A[0].shape[0]
    A2 = np.zeros(A.shape)
    frames=A.shape[0]
    print "Removing Drift"
    for i in range(A.shape[1]):
        A2[0][i][0] = A[0][i][0]
        A2[0][i][1] = A[0][i][1]
        A2[0][i][2] = A[0][i][2]
    for k in range(frames-1):
        for i in range(num):
            x = A[k+1][i][0] - Drift[k][0]
            y = A[k+1][i][1] - Drift[k][1]
            z = A[k+1][i][2] - Drift[k][2]
            x = boundary(x, L[0])
            y = boundary(y, L[1])
            z = boundary(z, L[2])
            A2[k+1][i][0] = x
            A2[k+1][i][1] = y
            A2[k+1][i][2] = z
    print "writing new file"
    f = open('drift.xyz','w')
    if index == False:
        for k in range(frames):
            f.write(('%i \nAtoms\n')%(A2.shape[1]))
            for j in range(A2[k].shape[0]):
                if j < A2[k].shape[0]/2:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',A2[k][j][0],A2[k][j][1],A2[k][j][2]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',A2[k][j][0],A2[k][j][1],A2[k][j][2]))
    else:
        for k in range(frames):
            f.write(('%i \nAtoms\n')%(A2.shape[1]))
            for j in range(A2[k].shape[0]):
                f.write(('%c %.2f %.2f %.2f \n')%(index[j],A2[k][j][0],A2[k][j][1],A2[k][j][2]))

    f.close()
    return A2
