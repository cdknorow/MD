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

## \brief Finds the spherical coordinates of a vector over timesteps 
# for a single Major particle under condtions where the box is changing size
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
def single_rotation_box(M,Z,L,index=0):
    theta = []
    phi = []
    sides = Z.shape[1]/M.shape[1]
    v = np.zeros((M.shape[0],3))
    for k in range(M.shape[0]):
        v[k] = points.vector1d(M[k][index],Z[k][index*sides],L[k])
    spher = xyz_to_spherical(v)
    theta = spher[:,1]
    phi = spher[:,2]
    return theta, phi

## \brief Finds the spherical coordinates of a vector over timesteps 
# for a single Major particle
#
# finds the theta coordinates where the center of the simulation box is
# the (0,0,0) axis. 
#
# \returns array theta coordinate at each time step
#
# \param M Matrix of Major particles M[frames][atoms][x,y,z]
# \param Z Matrix of sides of each central particle [frames][atoms][x,y,z] 
# \param sides number of sides for each Major particle
# \param index index of particle to look at
#
def rotation_angle(M,Z,L,index=0):
    theta = []
    v1 = []
    sides = Z.shape[1]/M.shape[1]
    for i in range(M.shape[0]):
        v1.append(points.vector1d(M[i][index],Z[i][index*sides],L[i]))
    for i in range(M.shape[0]-1):
        theta.append(points.angle_between(v1[i],v1[i+1]))
    return theta
## \brief Finds the spherical coordinates of a vector between t0 and t 
# for a single Major particle
#
# finds the theta coordinates where the center of the simulation box is
# the (0,0,0) axis. 
#
# \returns array theta coordinate at each time step
#
# \param M Matrix of Major particles M[frames][atoms][x,y,z]
# \param Z Matrix of sides of each central particle [frames][atoms][x,y,z] 
# \param sides number of sides for each Major particle
# \param index index of particle to look at
#
def diffusion(M,Z,L,index=0,t0=0,tfinal=1):
    theta = []
    v1 = []
    sides = Z.shape[1]/M.shape[1]
    for i in range(t0,tfinal):
        v1.append(points.unit(points.vector1d(M[i][index],Z[i][index*sides],L[i])))
    for i in range(tfinal-t0):
        theta.append(np.dot(v1[0],v1[i]))
    return theta
