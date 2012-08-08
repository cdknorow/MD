## \package MD.analysis.msd 
# \brief This module msd of atoms and takes into account drift of system by
# subtracting out center of mass.
#

import numpy as np
import MD.base.points as points


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
# \param A matrix of [frames][atoms][x,y,z] # \param L length of box
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
# \brief Find if the particle has crossed the boundary between timesteps
# and store that information.
#
# \returns updated information about crossing the boundary
#
# \param p1 point one
# \param p2 point two
# \param B stores information for how many times boundary has been crossed
#
def boundary_cross(p1,p2,L,B):
    d=p2-p1
    #if the particle has crossed the boundary
    #we add to 1 to B if it has crossed back we subtract from B
    if d[0] > L[0] / 2.0:
        B[0] += 1
    if d[0] <= -L[0] / 2.0:
        B[0] -= 1
    if d[1] > L[1] / 2.0:
        B[1] += 1
    if d[1] <= -L[1] / 2.0:
        B[1] -= 1
    if d[2] > L[2] / 2.0:
        B[2] += 1
    if d[2] <= -L[2] / 2.0:
        B[2] -= 1
    return B

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
# \param A matrix of [frames][atoms][x,y,z]  # \param L length of box
# \param step time step between dump files
#
def msd(A,L,step=1):
    #find the drift
    print A.shape
    Drift = drift_position(A,L)[0]
    #find msd
    num = A[0].shape[0]
    print num
    frames=A.shape[0]
    print "Calculating MSD"
    MSD=np.zeros((frames-1))
    B = np.zeros((num,3))
    for k in range(frames-1):
        r_sum=0
        for i in range(num):
            B[i] = boundary_cross(A[k][i], A[k+1][i], L, B[i])
            d = points.dist_np(A[0][i][0:3] + Drift[k], A[k+1][i][0:3] - B[i] * L)
            if d[0] > 1:
                r_sum = d[0]**2 / num + r_sum
        MSD[k] = r_sum
    x = np.arange(0, A.shape[0] * step - step, step)
    return x,MSD
## \brief Find the total displacemnt for a set of particles away from the center point
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
def diff(A,L,step=50000*5):
    #find the drift
    print A.shape
    Drift = drift_position(A,L)[0]
    #find msd
    num = A[0].shape[0]
    print num
    print num
    frames=A.shape[0]
    print "Calculating MSD"
    MSD=np.zeros((frames-1))
    r_sum = 0
    for k in range(0,frames-1):
        for i in range(num):
            d = points.dist(A[0][i][0:3], A[k][i][0:3], L)
            r_sum = abs(d[0])  + r_sum
        MSD[k] = r_sum
    x = np.arange(0, A.shape[0] * step - step, step)
    return x,MSD
## \brief Find the msd for a set of particles away from the center point
# with no drift
#
# msd is defined as 
# \f{eqnarray*}
# MSD(t) = \sum_{i=0}^N(r(0)_{i}-r(t)_i)za{2} 
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
def msd_no_drift(A,L,step=1):
    #find the drift
    #find msd
    num = A[0].shape[0]
    frames=A.shape[0]
    MSD=np.zeros((frames-1))
    B = np.zeros((num,3))
    for k in range(frames-1):
        r_sum=0
        for i in range(num):
            B[i] = boundary_cross(A[k][i], A[k+1][i], L, B[i])
            rt = A[k+1][i][0:3] - B[i] * L[0]
            ri = A[0][i][0:3]
            r = (ri-rt)
            r_sum += np.dot(r,r) / float(num)
        MSD[k] = r_sum**0.5
    x = np.arange(0, A.shape[0] * step - step, step)
    return x,MSD
