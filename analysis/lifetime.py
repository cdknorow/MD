## \package MD.analysis.lifetime
# \brief finds lifetime of connections made by groups of particles
#
# A connection is defined by two particles being within a distance rcut 
# simulation.
#

import numpy as np
from MD.base import points

## \internal
# \brief Check to see if particles are still connected.
#
# \returns particles that are still connected
# 
# \param A Matrix of Particles of group A [atoms][x,y,z]
# \param B Matrix of particles of group B [atoms][x,y,z]
# \param connected array of connected [index][con_A,con_B]
# \param L box length for periodic boundaries
# \param rcut distance with which particles are defined as connected
#
def still_connected(A,B,connected,L,rcut):
    con=[]
    for i in range(len(connected)):
        if points.dist(A[connected[i][0]],B[connected[i][1]],L)[0]<rcut:
            con.append(connected[i])
    return con


## \brief Find the lifetime of particle connections.
#
# \returns number of remaining connections after initial frame
# 
# \param A Matrix of Particles of group A [frames][atoms][x,y,z]
# \param B Matrix of particles of group B [frames][atoms][x,y,z]
# \param L box length for periodic boundaries
# \param rcut distance with which particles are defined as connected
# \param dh how many frames to skip between each check
def lifetime(A,B,L,rcut=1.0,dh=1):
    counter=[]
    #Find the first connection list
    master_connections=[]
    connected=[]
    for i in range(A.shape[1]):
        for j in range(B.shape[1]):
            if points.dist(A[0][i],B[0][j],L)[0]<rcut:
                master_connections.append([i,j])
    #Check to see if those connections still exist in each frame
    for k in range(0,A.shape[0]):
        #if they are connected there will be a 0 if not a 1
        still = still_connected(A[k],B[k],master_connections,L,rcut)
        master_connections=still
        connected.append(len(still))
        print connected[k/dh]
    connected = np.array(connected)/float(connected[0])
    return connected
