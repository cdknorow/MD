## \package MD.analysis.rotation 
# \brief This module is used to find rotations of particles in time
# 

import numpy as np
import MD.base.points as points
from MD.base.spherical import xyz_to_spherical

## \internal
# \brief Finds the spherical coordinates of a vector over
# timesteps
# 
# Note this method is not complete
#
# \returns array containing arrays of spherical coordanites
# at each timestep [frames][
# 
# \param M Matrix of Major particles M[frames][atoms][x,y,z]
# \param Z Matrix of sides of each central particle [frames][atoms][x,y,z] 
# \param sides number of sides for each Major particle
def rotations(M,Z,L,sides=6):
    theta = []
    phi = []
    zrange = range(0,Z.shape[0],sides)
    for i in range(M.shape[0]):
        v1 = points.vector(M[k],Z[k],L)
        spher = xyz_to_spherical(v1)
        theta.append(spher[1])
        phi.append(spher[2])
    return theta,phi

## \brief Finds the spherical coordinates of a vector over timesteps 
# for a single Major particle
#
# finds the theta and phi coordinates where the center of the simulation box is
# the (0,0,0) axis. 
#
# \returns array theta coordinate at each time step
# \returns array phi coordinate at each time step
#
# \param M Matrix of Major particles M[frames][atoms][x,y,z]
# \param Z Matrix of sides of each central particle [frames][atoms][x,y,z] 
# \param sides number of sides for each Major particle
# \param index index of particle to look at
#
def single_rotation(M,Z,L,index=0):
    theta = []
    phi = []
    sides = Z.shape[1]/M.shape[1]
    v1 = points.vector(M[:,index],Z[:,index*sides],L)
    spher = xyz_to_spherical(v1)
    theta = spher[:,1]
    phi = spher[:,2]
    return theta, phi

