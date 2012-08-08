## \package MD.base.points 
# \brief This module is methods for vector and point operations between
# particles. 
#
# This method also containes useful methods for finding the intersection, union,
# etc. It is used by several of the classes in MD and is meant to make code
# cleaner to read and get rid of repition. 
#
# most of the code is also writen to vectorize operations so that numpy arrays
# can be used to speed things up.


import numpy as np
from math import sqrt as sqrt
import copy
##\internal
#\brief union of elements in two lists (1d array)
# 
# \returns a combination of elements wihtout any repeats
# \param list1 1d array of elements
# \param list2 1d array of elements
def union(list1,list2):
    union=list1+filter(lambda x:x not in list1,list2)
    return union
##\internal
#\brief intersection of elements in two lists (1d array)
# 
# \returns elements which are both in list1 and list2
#
# \param list1 1d array of elements
# \param list2 1d array of elements
def intersection(list1,list2):
    intersection=filter(lambda x:x in list1,list2)
    return intersection
##\internal
#\brief difference of elements in two lists (1d array)
# 
# \returns elements which are in list1 but not list2
#
# \param list1 1d array of elements
# \param list2 1d array of elements
def difference(list1,list2):
    difference=filter(lambda x:x not in list2,list1)
    return difference
##\internal
#\brief distinct elements in two lists (1d array)
# 
# \returns elements which are in list1(2) but not list2(1)
#
# \param list1 1d array of elements
# \param list2 1d array of elements
def distinct(list1,list2):
    distinct=filter(lambda x:x not in list2,list1)+filter(lambda x:x not in
            list1,list2)
    return distinct
##\internal
#\brief unique elements in a single list(1d array)
# 
# \returns unique elements of list1
#
# \param list1 1d array of elements
def unique(list1):
    output = []
    for x in list1:
        if x not in output:
            output.append(x)
    return output
##\internal
# \brief looks through elements of array and returns a unique array among them

# \returns a unique array given many different arrays
# 
# \param array of arrays 
def unique_array(list1):
    for x in list1:
        x.sort()
    output=[]
    output.append(copy.deepcopy(list1[0]))
    for x in list1:
        count=0
        for j in range(len(output)):
            if x == output[j]:
                count+=1
        if count==0:
            output.append(x)
    return output
## \internal
# \brief gives a unit vector given 2 points
# 
# \returns unit vector
# \param p1 x,y,z of point 1
# \param p2 x,y,z of point 2
def unit_vect(p1,p2):
    d=p1-p2
    mag = sqrt(np.dot(d,d))
    unit = d/mag
    return unit
## \internal
# \brief gives a unit vector for a vector
# 
# \returns unit vector
# \param d x,y,z
def unit(d):
    mag = sqrt(np.dot(d,d))
    unit = d/mag
    return unit
## \internal
# \brief Applies boundary conditions for a point
#
# returns x corrected for boundary conditions
#
# \param x coordinate
# \param L box length of periodic boundaries
def apply_boundary(x,L):
    if x > L/2.0:
        x = x - L
    if x <- L / 2.0:
        x = x + L
    return x
## \internal
# \brief Applies boundary conditions for a vector
#
# returns x corrected for boundary conditions
#
# \param x coordinate
# \param L box length of periodic boundaries
def apply_vboundary(v,L):
    v[0] = apply_boundary(v[0],L)
    v[1] = apply_boundary(v[1],L)
    v[2] = apply_boundary(v[2],L)
    return v
## \internal 
# \brieffind the distances between two points whithout periodic boundary conditons
# returns distance between two points without using periodic boundaries
# 
# \param p1 x,y,z of points one
# \param p2 x,y,z of points two
def dist_np(p1,p2):
    #distance vector
    d=p2-p1
    #distance between p1,p2
    distance_scalar = np.sqrt(np.dot(d,d))
    return distance_scalar,d
## \internal
# \brief find the distances between two points while applying periodic boundary conditons
#
# \returns distance between two points or numpy arrays of p1,p2
# \returns vector between to points or numpy arrays of p1,p2
#
# \param p1 x,y,z of points one
# \param p2 x,y,z of points two
# \param L box length
def dist(p1,p2,L):
    #distance vector
    d=p2-p1
    #apply boundary conditions
    d[0]=apply_boundary(d[0],L[0])
    d[1]=apply_boundary(d[1],L[1])
    d[2]=apply_boundary(d[2],L[2])
    #distance scalar
    distance = sqrt(np.dot(d,d))
    return distance,d
## \internal
# \brief midpoint between two points where p1 is the origin 
#
# \returns vector to midpoint from origin
#
# \param p1 x,y,z of points one or numpy arrays
# \param p2 x,y,z of points two or numpy arrays
# \param L box length
def midpoint(p1,p2,L):
    d=(p2-p1)
    ##apply boundary conditions
    shape=d.shape
    d=map(lambda x: apply_boundary(x,L),d.ravel()) 
    d=np.array(d)
    d.resize(shape)
    d/=2
    return d 
## \internal
# \brief vector between two points including periodic boundaries 
#
# \returns vector to midpoint from origin
#
# \param p1 x,y,z of points one or numpy arrays
# \param p2 x,y,z of points two or numpy arrays
# \param L box length
def vector(p1,p2,L):
    v=(p2-p1)
    ##apply boundary conditions
    shape=v.shape
    print v.shape
    v.transpose()[0]=map(lambda x: apply_boundary(x,L[0]),v.transpose()[0])
    v.transpose()[1]=map(lambda x: apply_boundary(x,L[1]),v.transpose()[1])
    v.transpose()[2]=map(lambda x: apply_boundary(x,L[2]),v.transpose()[2])
    v=np.array(v)
    v.resize(shape)
    return v
## \internal
# \brief vector between two points including periodic boundaries 
#
# \returns vector to midpoint from origin
#
# \param p1 x,y,z of points one or numpy arrays
# \param p2 x,y,z of points two or numpy arrays
# \param L box length
def vector1d(p1,p2,L):
    v=(p2-p1)
    ##apply boundary conditions
    v[0]=apply_boundary(v[0],L[0])
    v[1]=apply_boundary(v[1],L[1])
    v[2]=apply_boundary(v[2],L[2])
    v=np.array(v)
    return v
## \internal
# \brief projection of vector B onto vector A 
#
# v = a.b/|a| * a/|a|
#
# \returns the projection of B onto A  
# \param A vector v1 
# \param B vector v2
def projection(A,B):
    dot = np.dot(A,B)
    unit = (np.dot(A,A))
    return A*dot/unit
## \internal
# \brief magnitude of vector A 
#
# v = a.b/|a| * a/|a|
#
# \returns the magnitude of vector A  
# \param A vector v1
def magnitude(A):
    return sqrt(np.dot(A,A))
## \internal
# \brief theta between two vectors v1,v2
#
# \returns theta between v1 and v2 (radians)  
# \param v1 vector v1 
# \param v2 vector v2
def angle_between(v1,v2):
    return np.arccos(np.dot(v1,v2)/(magnitude(v1)*magnitude(v2)))
## \internal
# \breif projection of vector onto a plane
#
# return vector projected onto a plane
#
# \param v1 vector to be projected with origin at 0,0,0
# \param p1 vector of plane with origin at 0,0,0 
def vector_plane_projection(v,p):
    return np.cross(p,np.cross(v,p)/magnitude(p))/magnitude(p)
