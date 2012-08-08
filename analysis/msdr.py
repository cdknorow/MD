## \package MD.analysis.msdr
# \brief Find the deviation of every particle from its point on 
# the recipricol lattice.
#

import os
import math
import numpy as np
import MD.base.points as points

## \brief given the reconstructed lattice and the x,y,z of the real coordinates
# find the mean square displacement from the reconstructed and the real lattice
#
# \return msdr 
# 
# \param A Major points in a frame A[atoms][ x, y, z]
# \param lattice Major points in a frame lattice[atoms][ x, y, z] of
# reconstructed lattice
# \param L box length
def msd_r(A,lattice,L):
    #First find the closest to the real coordinate
    msd=np.zeros((A.shape[0],2))
    msd+=L[0]+10
    for i in range(A.shape[0]):
        for j in range(lattice.shape[0]):
            D=points.dist(A[i],lattice[j],L)[0]
            #if the distance is shorter store the point and the distance
            if abs(msd[i][1])>abs(D):
                msd[i][0]=j
                msd[i][1]=D
    #find the total msd
    MSD=0
    for i in range(msd.shape[0]):
        MSD+=msd[i][1]**2
    MSD_total= MSD/A.shape[0]
    print 'msd for choice',MSD_total
    return MSD_total

