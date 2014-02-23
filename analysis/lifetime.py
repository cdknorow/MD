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
def still_connected(A,B):
    con=[[],[]]
    for i in range(len(A[0])):
        try:
            index = B[0][0].index(A[0][i])
            if A[1][i] == B[0][1][index]:
                con[0].append(A[0][i])
                con[1].append(A[1][i])
            else:
                try:
                    index2 = B[1][0].index(A[0][i])
                    if A[1][i] == B[1][1][index2]:
                        con[0].append(A[0][i])
                        con[1].append(A[1][i])
                        print index2
                except:
                    pass
        except:
            pass
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
def lifetime(connections):
    counter=[]
    #Find the first connection list
    life=[]
    master_connections = connections[0]
    #Check to see if those connections still exist in each frame
    for k in range(len(connections)):
        #if they are connected there will be a 0 if not a 1
        master_connections =  still_connected(master_connections,connections[k:k+1])
        life.append(len(master_connections[0]))
    return life
