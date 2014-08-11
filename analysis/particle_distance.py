## \package MD.analysis.particle_distance 
# \brief This module is used to find the distance between many particles
#
# This is mainly used for generating a distribution of distances between
# particles

import numpy as np

# \brief find the distances between Points in Two arrarys arrays must have same
# elements 
# 
# \returns array of distances
#
# \param A matrix of A points for a A[atoms][x,y,z]
# \param B matrix of B points for a B[atoms][x,y,z]
# \param L lenght of box 
def particle_distance(A,B,L):
    d = np.zeros(A.shape[0]*A.shape[1]*A.shape[1])
    counter = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            r = np.abs(B[i] - A[i][j])
            r = np.where(r > L[0]/2, L[0] - r, r)
            d[counter:counter+A.shape[1]] = np.sqrt((r**2).sum(axis=-1))
            counter += A.shape[1]
    return d

#def particle_distance(A,B,L,rmax=1000):
#    d = []
#    for i in range(A.shape[0]):
#        for j in range(A.shape[1]):
#            for k in range(B.shape[1]):
#                r = dist(A[i][j],B[i][k],L)[0]
#                if r > 0 and r < rmax:
#                    d.append(r)
#    return d

# \brief attempts to find the first peak of a distance distribution
#
# \returns first peak it finds 
#
# \param A matrix of A points for a A[atoms][x,y,z]
# \param B matrix of B points for a B[atoms][x,y,z]
# \param L lenght of box 
# \param cut cutoff for declaring the first peak
# \param filter filter to not call it a peak below certain value
def min_particle_distance(A,B,L,cut=10,verbose=False):
    d = particle_distance(A,B,L)
    hist,x = np.histogram(d,bins=75)
    #look a peak that is always greater than the next point 
    #for at least 5 consecutivej points
    hist = hist[1:]
    x = x[1:]
    first_peak = 0
    count = 0
    for i in range(x.shape[0]):
        if hist[first_peak] > hist[i]:
            count+=1
        else:
            first_peak = i
            count = 0
        if count > cut:
            if verbose:
                print 'first peak found',x[first_peak]
                print hist
                print x
            break
    return x[first_peak]

# \brief find the smallest distances between a point and points in an array array
# elements 
# 
# \returns array of distances
#
# \param P points for [x,y,z]
# \param B matrix of B points for a B[atoms][x,y,z]
# \param L lenght of box 
def min_point_distance(P,B,L):
    small = 100
    index = 0
    for i in range(B.shape[0]):
        d = dist(P,B[i],L)
        if d[0] < small:
            small = d[0]
            vector = d[1] 
            index = i
    return small,vector,index


if __name__ == '__main__':
    A = np.zeros((1,5,3))
    A[0][1] += 1 
    A[0][2] += 1 
    A[0][3] += -5 
    A[0][4] += 5 
    print A
    print particle_distance(A, A,[10,10,10])
