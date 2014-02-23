import os
import random
import numpy as np
import points

def sort_clusters(CM):
    clusters = []
    c = np.array([0,0,0])
    c_rank = []
    for i in range(len(CM)):
        d = points.dist_np(c,CM[i])[0]
        c_rank.append([d,i])
    c_rank.sort()
    for i in c_rank:
        clusters.append(CM[i[1]])
    return clusters

#\brief returns neighbor indexes
#
#\param N array containing x,y,z coordinates
#\param i index of particle
#\param rcut cutoff distance for finding a neighbor
#
#\returns array containing index of neighboring particles
def neighbors(N,i,rcut):
    neighbors = []
    for index in range(len(N)):
        if i !=index:
            if points.dist_np(N[index],N[i])[0] < rcut:
                neighbors.append(index)
    return neighbors

#\brief returns center of mass of clusters (non periodic)
#
#\param N array containing x,y,z coordinates
#\param Cluster arrays of clustered particles
#
#\returns x,y,z of cluster coordinates for center of mass of cluster
def average_clusters(N,Cluster):
    CM = []
    for i in range(len(Cluster)):
        v_sum = np.zeros((3))
        count = 0
        for v in Cluster[i]:
            v_sum += N[v]
            count += 1
        CM.append(v_sum/count)
    return sort_clusters(CM)

#\brief returns center of mass of clusters (non periodic)
#
#\param N array containing x,y,z coordinates
#\param rcut cutoff distance for finding a neighbor
#
#\returns x,y,z of cluster coordinates
def cluster(N,rcut=1.5,cluster_cut=6):
    NN=[]
    for i in range(len(N)):
        NN.append(neighbors(N,i,rcut))
    for i in range(200):
        count = random.randint(0,len(NN)-1)
        N_new = []
        for i in range(len(NN)):
            if len(points.intersection(NN[count],NN[i]))>0:
                NN[count]=points.union(NN[count],NN[i])
            elif len(NN[i])>0:
                N_new.append(NN[i])
        N_new.append(NN[count])
        NN = N_new
    N_new = []
    for i in NN:
        if len(i) >cluster_cut:
            N_new.append(i)
    NN = N_new
    return average_clusters(N,NN)







