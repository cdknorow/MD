
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
cimport numpy as np
import math
import copy

# Fix Data Type for array
DTYPE = np.float
ITYPE = np.int
# Assign compile time Data Type for array
ctypedef np.float_t DTYPE_t
ctypedef np.int_t ITYPE_t

## \internal
# \brief Applies boundary conditions for a point
#
# returns x corrected for boundary conditions
#
# \param x coordinate
# \param L box length of periodic boundaries
def apply_boundary(float x, float L):
    if x > L:
        x = x - L
    if x <- L:
        x = x + L
    return x
## \internal
# \brief find the distances between two points while applying periodic boundary conditons
#
# \returns distance between two points or numpy arrays of p1,p2
# \returns vector between to points or numpy arrays of p1,p2
#
# \param p1 x,y,z of points one
# \param p2 x,y,z of points two
# \param L box length
def dist(np.ndarray[DTYPE_t, ndim=1] p1,
        np.ndarray[DTYPE_t, ndim=1] p2,
        np.ndarray[DTYPE_t, ndim=1] L):
    #distance vector
    cdef np.ndarray[DTYPE_t, ndim=1] d = np.zeros((3), dtype = DTYPE)
    cdef float distance
    d=p2-p1
    #apply boundary conditions
    d[0]=apply_boundary(d[0],L[0]/2.)
    d[1]=apply_boundary(d[1],L[1]/2.)
    d[2]=apply_boundary(d[2],L[2]/2.)
    return d
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
def drift_position(np.ndarray[DTYPE_t, ndim=3] A,
                   np.ndarray[DTYPE_t, ndim=1] L,
                   int start):
    cdef int num = A[0].shape[0]-start
    print "Calculating drift of position"
    cdef np.ndarray[DTYPE_t, ndim=2] D = np.zeros(((A.shape[0]-1),3), dtype = DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] dx = np.zeros((3), dtype = DTYPE)
    cdef float x_whole = 0
    cdef float y_whole = 0
    cdef float z_whole = 0
    cdef float x_sum, y_sum, z_sum
    cdef int i,k
    cdef int finish  = A.shape[0] - 1
    for k in range(start, finish):
        x_sum = 0; y_sum = 0; z_sum = 0
        for i in range(num):
            dx = dist(A[k][i], A[k+1][i], L)
            x_sum += dx[0]
            y_sum += dx[1]
            z_sum += dx[2]
        D[k][0] = (x_sum / num + x_whole)
        D[k][1] = (y_sum / num + y_whole)
        D[k][2] = (z_sum / num + z_whole)
        x_whole = x_sum / num + x_whole
        y_whole = y_sum / num + y_whole
        z_whole = z_sum / num + z_whole
    return  D

## \internal # \brief Find if the particle has crossed the boundary between timesteps
# and store that information.
#
# \returns updated information about crossing the boundary
#
# \param p1 point one
# \param p2 point two
# \param B stores information for how many times boundary has been crossed
#
def boundary_cross( np.ndarray[DTYPE_t, ndim=1] p1,
                    np.ndarray[DTYPE_t, ndim=1] p2,
                    np.ndarray[DTYPE_t, ndim=1] L,
                    np.ndarray[ITYPE_t, ndim=1] B):
    cdef np.ndarray[DTYPE_t, ndim=1] d = np.zeros((3), dtype = DTYPE)
    d=p2-p1
    #if the particle has crossed the boundary
    #we add to 1 to B if it has crossed back we subtract from B
    if d[0] > L[0]:
        B[0] += 1
    if d[0] <= -L[0]:
        B[0] -= 1
    if d[1] > L[1]:
        B[1] += 1
    if d[1] <= -L[1]:
        B[1] -= 1
    if d[2] > L[2]:
        B[2] += 1
    if d[2] <= -L[2]:
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
def int_scatter(np.ndarray[DTYPE_t, ndim=3] A,
                np.ndarray[DTYPE_t, ndim=1] L,
                np.ndarray[DTYPE_t, ndim=1] q,
                np.ndarray[DTYPE_t, ndim=1] x,
                int step=1):
    #find the drift
    cdef np.ndarray[DTYPE_t, ndim=2] Drift
    Drift = drift_position(A,L,x[0])
    #find msd
    cdef int N = A.shape[1]
    cdef int frames = A.shape[0]
    cdef int x_len = len(x) 
    print "Calculating Self-Intermediate Scattering Function"
    cdef np.ndarray[DTYPE_t, ndim=1] F = np.zeros((x_len), dtype = DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] B = np.zeros((N,3), dtype = ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] r = np.zeros((3), dtype = DTYPE)
    cdef int i, j, k
    for j,k in enumerate(x):
        for i in range(N):
            r = (A[x[0]][i][0:3] + Drift[k-1]) - (A[k][i][0:3] - B[i] * L)
            B[i] = boundary_cross(A[k][i], A[k+1][i], L/2.0, B[i])
            F[j] += math.cos(np.dot(q,-r))
        F[j] /= N
    return F
