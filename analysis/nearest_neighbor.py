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
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param index index of particle to find nearest neighbors of
# \param L length of box correspodning to periodic boundaries
# \param rcut cutoff distance for nearest neighbor search
#
# To calculated the nearest neighbors M, index, and L must be specified
# rcut will be set to a default.
#
#
def nearest_neighbors_index(M,index,L,rcut=20):
    neighbors = []
    vectors = []
    for i in range(M.shape[0]):
        if dist(M[i],M[index],L)[0] <= rcut and i!=index:
            neighbors.append(i)
    vectors = vector(M[index],M[neighbors],L)
    if neighbors == []:
        print "error No neighbors were found try reducing rcut"
    return neighbors, vectors

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
    neighbors = []
    vectors = []
    for i in range(M.shape[0]):
        if (dist(M[i],point,L)[0] <= rcut
        and np.array_equal(point,M[i])==False):
            neighbors.append(i)
    vectors = vector(point,M[neighbors],L)
    if neighbors == []:
        print "error No neighbors were found try increasing"
    return neighbors, vectors

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
## \brief Calculates the nearest neighbors to single particle.
# 
# \returns index of nearest neighbors
# \returns vectors to nearest neighbors
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param index index of particle to find nearest neighbors of
# \param L length of box correspodning to periodic boundaries
# \param rcut cutoff distance for nearest neighbor search
#
# To calculated the nearest neighbors M, index, and L must be specified
# rcut will be set to a default.
#
#
def second_nearest_neighbors_index(M,index,L,inrcut=20, outrcut=40):
    neighbors = []
    vectors = []
    for i in range(M.shape[0]):
        d = dist(M[i],M[index],L)[0] 
        if (d <= outrcut and  i != index and d <= inrcut):
            neighbors.append(i)
    vectors = vector(M[index],M[neighbors],L)
    if neighbors == []:
        print "error No neighbors were found try reducing rcut"
    return neighbors, vectors
## \brief finds the closest neighbors to a point x,y,z
# 
# \returns index of nearest neighbor
# \returns distance to nearest neighbor
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param point [x,y,z] point to find neihbors of
# \param L length of box correspodning to periodic boundaries
#
def ws_neighbors_point(M,point,L,index,rmin=6,rmax=10):
    neighbors_in = []
    neighbors_out = []
    for i in range(M.shape[0]):
        if dist(M[i],point,L)[0] <= rmin:
            neighbors_in.append(i)
        if dist(M[i],point,L)[0] <= rmax:
            neighbors_out.append(i)
    #if len(neighbors_in) != len(neighbors_out):
    #    print 'neighbors not equal', index
    #    print neighbors_in
    #    print neighbors_out
    #if len(neighbors_out) == 2:
    #    print 'possible intersticial', index
    return neighbors_in, neighbors_out

## \brief finds num neighbors to a point x,y,z within rcut
# 
# \returns index of nearest neighbor
# \returns distance to nearest neighbor
#
# \param M Numpy matrix of [frames][particles][x,y,z]
# \param point [x,y,z] point to find neihbors of
# \param L length of box correspodning to periodic boundaries
#
def rcut_neighbors_point(M,point,L,rcut=10):
    neighbors = []
    for i in range(M.shape[0]):
        p = dist(M[i],point,L)[0]
        if dist(M[i],point,L)[0] <= rcut :
            neighbors.append(i)
    if abs(point[0]) < 15 and abs(point[1]) < 15 and abs(point[2]) < 15:
        if len(neighbors) > 0:
            return len(neighbors)
        else:
            return 0
    if len(neighbors) > 0:
        return len(neighbors)
    return 0

