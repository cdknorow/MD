
# \package MD.analysis.int_scattering
# \brief This module finds the self intermediate scattering function peaks for a group of
# a trajectory at equillibrium.
#
# \b Example
# \code
#    A = Average(VW,L,n_start,n_frames)
#    stmp,qx = sf.sfactor(np.array([A]),L=L[n_start],l=l)
#    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
#    k = primitive_vectors[S.index(max(S))]
#    x,F = int_scatter(VW,L,k,x)
# \endcode
#
# S is the S(q), Q is mag(q_vec), and primitive vectors are q vectors

import os
import numpy as np
import math
import copy
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
def drift_position(A,L,start):
    num = A[0].shape[0]-start
    print "Calculating drift of position"
    D = np.zeros(((A.shape[0]-1),3))
    Dsum = np.zeros(((A.shape[0]-1),1))
    x_whole = 0; y_whole = 0; z_whole = 0
    for k in range(start,A.shape[0] - 1):
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

## \internal # \brief Find if the particle has crossed the boundary between timesteps
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

## \brief Find the self intermediate scattering function F_s(k,t)
#
# F_s(k,t) is defined as 
# \f{eqnarray*}
# F_s(k,t) = 1/N\sum_{i=0}^N(<exp[ k dot (r(0)_{i}-r(t)_i)]>} 
# \f}
#where k is the first peak of the structure facture S(q)
# r vector of the particle and ri is ith particle in frame k. 
#
#the drift of the system is also subtracted away to get an accurate measure 
# of the movement of the particles
# 
# \returns time steps and corresponding F(k,t)
#
# \param A matrix of [frames][atoms][x,y,z]  # \param L length of box
# \param step time step between dump files
# \k S(k) peak from structure factor
#
def int_scatter(A,L,q,x,step=1):
    #find the drift
    Drift = drift_position(A,L,x[0])[0]
    #find msd
    N = A.shape[1]
    frames = A.shape[0]
    print "Calculating Self-Intermediate Scattering Function"
    F=np.zeros((len(x)))
    B = np.zeros((N,3))
    for j,k in enumerate(x):
        for i in range(N):
            r = (A[x[0]][i][0:3] + Drift[k-1]) - (A[k][i][0:3] - B[i] * L)
            B[i] = boundary_cross(A[k][i], A[k+1][i], L, B[i])
            F[j] += math.cos(np.dot(q,-r))
        F[j] /= N
    return F
