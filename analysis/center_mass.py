## \package MD.analysis.center_mass 
# \brief This module is used to find the center of mass of particles 
#
# A good use of this is to create major particles from polymers in order to
# better make use of the anlysis tools in this package
#
# Note this method has not been fully implemented

import numpy as np
import MD.base.points as points

## \internal
# \brief apply a boundary condition on particle distances 
#
# We apply a boundary condition at L/4 to bring everything into the center of
# the box
#
# \returns x coordinate after applying boundary conditions(if applied)
# \returns applied flag(set to true if boundary conditions were applied)
#
# \param x coordinate to analyze
# \param L Box length
##########################################
def cm_boundary(x,L,applied = False):
    if x>L/4.0:
        x=x-L/2.0
        applied=True
    if x<-L/4.0:
        x=x+L/2.0
        applied=True
    return x,applied

##\internal
# \brief correct for the cm_boundary
def reset_boundary(x,L,applied):
    if x>0:
        x=x-L/2.0
    if x<0:
        x=x+L/2.0
    return x
## \brief Finds the CM of a cluster of particles given 
#  
# \returns np.array containing xyz of cm for each cluster
#
# \param M Matrix of points at a given frame
# \param cluster array of arrays of cluster index points       
# \param L box size
# 
##########################################
def CM(M,cluster,L):
    cm = np.zeros((len(cluster),3))
    for i in range(len(cluster)):
        #This will let us know if we used the cm_bounday
        #so that we can correct the Cm back to its proper value
        # after the calculation
        x_set = False
        y_set = False
        z_set = False
        #Find CM 
        for j in range(len(cluster[i])):
            xcm, x_set =  cm_boundary(M[cluster[i]][0], L, x_set)[0]:
            ycm, y_set =  cm_boundary(M[cluster[i]][1], L, y_set)[0]:
            zcm, z_set =  cm_boundary(M[cluster[i]][2], L, z_set)[0]:
            x_total += xcm
            y_total += ycm
            z_total += zcm
        if x_set:
            x_total = rest_boundary(x_total)
        if y_set:
            y_total = rest_boundary(y_total)
        if z_set:
            z_total = rest_boundary(z_total)
        cm[i][0] = x_total
        cm[i][1] = y_total
        cm[i][2] = z_total
    return cm
