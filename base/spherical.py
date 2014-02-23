## \package MD.base.spherical 
# \brief This module is methods for dealing with spherical math 
#
#
# (note) most of the code is also writen to vectorize operations so that numpy arrays
# can be used to speed things up.

import numpy as np
import math

##############################
#takes xyz points and makes them into spherical coordinates
#thea => 0,2pi
#phi => 0,pi
############################
def xyz_to_spherical(xyz):
        spher=np.zeros(xyz.shape)
        xy = xyz[:,0]**2 + xyz[:,1]**2
        #r,theta,phi
        spher[:,0] = np.sqrt(xy + xyz[:,2]**2)
        #theta (azimuth)
        spher[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
        #phi (inclination/polar)
        spher[:,2] = np.arccos(xyz[:,2]/spher[:,0])
        for i in range(spher.shape[0]):
            if spher[i,1]<0:
                spher[i,1]+=2*math.pi
        return spher
##########################################
#spherical harmoincs using scipy
#M is array r,theta,phi
#n,m are l amd m
##########################################
def spherical_harmonics(M,l,m):
    import scipy.special as special
    spher=special.sph_harm(m,l,M[:,1],M[:,2])
    return spher
#######################
#complex vector for particle i
#M is array of spherical cordinates from particle to neighbors
#####################
def qlm(M,l,m):
    qlm=(spherical_harmonics(M,l,m)).sum()
    return np.round(qlm/M.shape[0],decimals=8)
################################
#Ql 
#returns ql for l 
#################################
def ql(qlm,l):
    #other stuff
    ql =(4*math.pi/(2*l+1))*((qlm*qlm.conjugate()).sum())
    #ql =((qlm*qlm.conjugate()).sum())
    return ql**0.5
################################
#Sij sum over m for l in qlm 
#two particles are connected if Sij>0.5
#return Sij value between q1 particles
#and each of its neighbors as array
#################################
def Sij(q1,q_neighbors):
    #Qlm[i]=np.array(ql-m:ql+m) for particle i
    #Sij[i][:]=Sij value between x and particle j
    S=np.zeros((q_neighbors.shape[0],1),dtype=complex)
    for i in range(S.shape[0]):
        #S[i]=((q1*q_neighbors[i].conjugate()+q1.conjugate()*q_neighbors[i])/2.0).sum()
        S[i]=(q1*q_neighbors[i].conjugate()).sum()
    return S


if __name__ == "__main__":
    M = np.array([[0,0,1],[0,0,-1],[0,-1,0],[1,1,1],[1,1,0],[0,1,1]])
    spher = xyz_to_spherical(M)
    print spher
    print np.round(spherical_harmonics(spher,1,1),decimals=6)

