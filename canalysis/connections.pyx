## \package MD.analysis.connections 
# \brief This module is used to find connections between particle over a series
# of frames
# 
# Methods for finding connections between particles with one or more sites that
# can make connections with other particles in the system. A connecing site will
# be defined as a place where two particles come withing a certain cutoff
# distance of each other. 
# 

import numpy as np
cimport numpy as np
cimport cython


FTYPE = np.float
ctypedef np.float_t FTYPE_t


## \brief apply boundary conditions to a vector
# \param p1, p2 two points  
# \param L Length of Box 
def dist(np.ndarray[FTYPE_t, ndim=1] p1, np.ndarray[FTYPE_t, ndim=1] p2, np.ndarray[FTYPE_t, ndim=1] L):
    cdef float d, x, y, z
    #distance vector
    p1 = p2 - p1
    #apply boundary conditions
    x = p1[0]
    if x>L[0]/2.0:
        x=x-L[0]
    if x<-L[0]/2.0:
        x=x+L[0]
    y = p1[1]
    if y>L[1]/2.0:
        y=y-L[1]
    if y<-L[1]/2.0:
        y=y+L[1]
    z = p1[2]
    if z>L[2]/2.0:
        z=z-L[2]
    if z<-L[2]/2.0:
        z=z+L[2]
    p1[0] = x
    p1[1] = y
    p1[2] = z
    #distance scalar
    d = np.sqrt(np.dot(p1,p1))
    return d


## \brief Find connections within rcut of each other and return connection numbers
# \param A numpy matrix of [frames][atoms][x, y, z]
# \param B numpy matrix of [frames][atoms][x, y, z]
# \param L Length of Box 
# \param rcut cutoff distance for a connectio
def connections(np.ndarray[FTYPE_t, ndim=3] A, np.ndarray[FTYPE_t, ndim=3] B,
        np.ndarray[FTYPE_t, ndim=1] L, float rcut = 1.5):
    counter=[]
    cdef int frames = A.shape[0]
    cdef int amax = A.shape[1]
    cdef int bmax = B.shape[1]
    cdef float x, y, z, d
    cdef np.ndarray[FTYPE_t, ndim=1] p = np.zeros([3],dtype=np.float)
    for k in range(frames):
        print k
        count = []
        for i in range(amax):
            for j in range(bmax):
                d = dist(A[k][i],B[k][j],L)
                if d < rcut:
                    count.append([i,j+amax])
        counter.append(count)
    return counter

