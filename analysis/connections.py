## \package MD.analysis.connections 
# \brief This module is used to find connections between particle over a series
# of frames
# 
# Methods for finding connections between particles with one or more sites that
# can make connections with other particles in the system. A connecing site will
# be defined as a place where two particles come withing a certain cutoff
# distance of each other. 
# 

from MD.base import points
from MD.base.points import dist
import numpy as np

## \brief Find connections within rcut of each other and return connection numbers
# \param A numpy matrix of [frames][atoms][x, y, z]
# \param B numpy matrix of [frames][atoms][x, y, z]
# \param L Length of Box 
# \param rcut cutoff distance for a connection
def connections(A,B,L,rcut=1.5):
    counter=[]
    for k in range(A.shape[0]):
        count = []
        print 'step',k
        for i in range(A.shape[1]):
            for j in range(B.shape[1]):
                d = dist(A[k][i],B[k][j],L)[0]
                if d<rcut:
                    count.append([i,j+A.shape[1]])
        counter.append(count)
    return counter

## \brief Find connections within rcut of each other and return connection numbers
# \param A numpy matrix of [frames][atoms][x, y, z]
# \param B numpy matrix of [frames][atoms][x, y, z]
# \param L Length of Box 
# \param rcut cutoff distance for a connection
def connections_min_max(A,B,L,rcut=1.5,rmin=1):
    counter=[]
    for k in range(A.shape[0]):
        count = []
        print 'step',k
        for i in range(A.shape[1]):
            for j in range(B.shape[1]):
                d = dist(A[k][i],B[k][j],L)[0]
                if (d<rcut) & (d>rmin):
                    count.append([i,j+A.shape[1]])
        counter.append(count)
    return counter
## \brief Find the number of connections to each particle
# if you have system whcih has particles which can have many connections, then
# this will return the average number of connections per particle.
# 
# \returns array of connections at each frame [frame][number of connections]
#
# \param connections list of connections at each time step [frame][[p1,p6],[p3,p10]]
#   where p1 and p6 would be two particles with a connection
# \param N  the number of of central points which can be connected too ie.
# spheres with dna attatched where the dna hybridize to other dna N would be the
# number of spheres
#########################################
def num_connections(connections,N):
    con=[]
    print 'finding # of connections'
    for i in range(len(connections)):
        con.append(len(connections[i][0])/float(N)*2)
    return np.array(con)
## \brief find number of connections where a linker is involved
# This is a case where there is a linking connection which can link two
# particles in some cases we would like to know if the linker is connected at
# on both/one/none of its connecting sites.
#
# \returns number of linkers with two sites connected
# \returns number of linkers with a single site connected
# \returns list of linkers with a single site connected
# \returns list of linkers with a both site connected
#
# \param con_A = list of connections at each time step for site A of the linker 
# \param con_B = list of connections at each time step for site B of the linker 
#########################################
def both_connections(con_A,con_B):
    l1 = []
    l2 = []
    #get the values for linekrs which have connectison
    #in an agreeable format
    for k in range(len(con_A)):
        a = []
        b = []
        for i in range(len(con_A[k])):
            a.append(con_A[k][i][1])
        for i in range(len(con_B[k])):
            b.append(con_B[k][i][1])
        l1.append(a)
        l2.append(b)
    #count the number that are unique and the same
    one=[]
    two=[]
    two_list=[]
    one_list=[]
    for k in range(len(l1)):
        two.append(len(points.intersection(l1[k],l2[k])))
        one.append(len(points.difference(l1[k],l2[k])))
    for k in range(len(l1)):
        two_list.append((points.intersection(l1[k],l2[k])))
        one_list.append((points.difference(l1[k],l2[k])))
    return one,two,one_list, two_list

## \brief find number of connections where a linker is a polymer
# This is a case where there is a linking connection which can link two
# particles in some cases we would like to know if the linker is connected at
# on both/one/none of its connecting sites.
#
# \returns number of linkers with two sites connected
# \returns number of linkers with a single site connected
# \returns list of linkers with a single site connected
# \returns list of linkers with a both site connected
#
# \param con_A = list of connections at each time step for site A of the linker 
#########################################
def both_connections_polymer(con_A,num):
    l1 = []
    l2 = []
    #get the values for linekrs which have connectison
    #in an agreeable format
    for k in range(len(con_A)):
        a = []
        b = []
        for i in range(len(con_A[k])):
            if con_A[k][i][1]%2 == 0:
                a.append(con_A[k][i][1]-num)
            if con_A[k][i][1]%2 == 1:
                b.append(con_A[k][i][1]-num)
        l1.append(a)
        l2.append(b)
    #count the number that are unique and the same
    one=[]
    two=[]
    two_list=[]
    one_list=[]
    for k in range(len(l1)):
        two.append(len(points.intersection(l1[k],l2[k])))
        one.append(len(points.difference(l1[k],l2[k])))
    for k in range(len(l1)):
        two_list.append((points.intersection(l1[k],l2[k])))
        one_list.append((points.difference(l1[k],l2[k])))
    return one,two,one_list, two_list

