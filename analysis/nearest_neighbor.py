## \package MD.analysis.nearest_neighbor
# \brief methods for finding nearest neighbors of particles in a
# simulation.
#
# This method containts two functions which can find the nearest 
# neighbors of a point or the nearest neighbor of a particle given
# its index in the matrix M.
#

import numpy as np
import MD.base.points as points
from MD.base.points import dist
from MD.base.points import vector


## \brief Calculates the nearest neighbors to single particle.
# 
# \returns index of nearest neighbors
# \returns vectors to nearest neighbors
#
# \param M Numpy matrix of [particles][x,y,z]
# \param index index of particle to find nearest neighbors of
# \param L length of box correspodning to periodic boundaries
# \param rcut cutoff distance for nearest neighbor search
#
# To calculated the nearest neighbors M, index, and L must be specified
# rcut will be set to a default.
#
#
def nearest_neighbors_index(M,index,L,rcut=20, vectors = False):
    if vectors:
        r = M - M[index]
        r = np.where(r > L[0]/2, r - L[0], r)
        r = np.where(r < -L[0]/2, r + L[0], r)
        d = np.sqrt((r**2).sum(axis=-1))
        return np.argwhere(d < rcut)[1:].tolist(), r[np.argwhere(d < rcut)[1:]]
    else:
        r = np.abs(M - M[index])
        r = np.where(r > L[0]/2, L[0] - r, r)
        d = np.sqrt((r**2).sum(axis=-1))
        return np.argwhere(d < rcut)[1:].tolist()


## \brief Finds x number of neighbors to single particle.
# 
# \returns index of nearest neighbors
# \returns vectors to nearest neighbors
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param index index of particle to find nearest neighbors of
# \param L length of box correspodning to periodic boundaries
# \param count number of neighbors to find
#
# To calculated the nearest neighbors M, index, and L must be specified
# rcut will be set to a default.
#

def count_neighbors_index(M,index,L,count=8,rcut=20, vectors = False):
    if vectors:
        r = M - M[index]
        r = np.where(r > L[0]/2, r - L[0] , r)
        r = np.where(r < -L[0]/2, r + L[0], r)
        d = np.sqrt((r**2).sum(axis=-1))
        b = filter(lambda x: d[x] < rcut,d.argsort()[1:count+1] ) 
        return b , r[b]
    else:
        r = np.abs(M - M[index])
        r = np.where(r > L[0]/2, L[0] - r, r)
        d = np.sqrt((r**2).sum(axis=-1))
        b = filter(lambda x: d[x] < rcut,d.argsort()[1:count+1] )
        return b


## \brief finds the nearest neighbors to a point x,y,z
# 
# \returns index of nearest neighbors
# \returns vectors to nearest neighbors
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param point [x,y,z] point to find neihbors of
# \param L length of box correspodning to periodic boundaries
# \param rcut cutoff distance for nearest neighbor search
#
# To calculated the nearest neighbors M, point, and L must be specified
# rcut will be set to a default.
#
#
def nearest_neighbors_point(M,point,L,rcut=20):
    r = np.abs(M - point)
    r = np.where(r > L[0]/2, L[0] - r, r)
    d = np.sqrt((r**2).sum(axis=-1))
    if vectors:
        return np.argwhere(d < rcut)[1:], r[np.argwhere(d < rcut)[1:]]
    return np.argwhere(d < rcut)[1:]


## \brief finds the closest neighbor to a point x,y,z
# 
# \returns index of nearest neighbor
# \returns distance to nearest neighbor
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param point [x,y,z] point to find neihbors of
# \param L length of box correspodning to periodic boundaries
#
def close_neighbors_point(M,point,L):
    neighbors = False
    lowest = 100
    for i in range(M.shape[0]):
        if dist(M[i],point,L)[0] != 0:
            D = dist(M[i],point,L)[0]
            if lowest > D:
                lowest = D
                neighbor = i
    return neighbor, lowest

## \brief finds the second closest neighbor to a point x,y,z
# 
# \returns index of nearest neighbor
# \returns distance to nearest neighbor
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param point [x,y,z] point to find neihbors of
# \param L length of box correspodning to periodic boundaries
#
def close_2neighbors_point(M,point,L):
    neighbors = []
    vectors = []
    lowest = 100
    neighbor2 = 0
    for i in range(M.shape[0]):
        if dist(M[i],point,L)[0] != 0:
            D = dist(M[i],point,L)[0]
            if lowest > D:
                lowest = D
                neighbor = i
                neighbor2 = neighbor
    try:
        return neighbor, neighbor2, lowest
    except:
        print "error No neighbors were found try increasing"

## \brief finds the closest point between two arrays
# 
# \returns index of closet points in array one
# \returns index of closet point in array two
#
# \param A1,A2 arrays of points you are interested 
# \param L length of box correspodning to periodic boundaries
#
def min_distance(A1,A2,L):
    num=0
    side1=0
    side2=0
    smallest=50
    for i in range(A1.shape[0]):
        for j in range(A2.shape[0]):
            num = dist(A1[i],A2[j],L)[0]
            if smallest > num:
                smallest = num
                side1 = i
                side2 = j
    return side1,side2

