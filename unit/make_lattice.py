## \package MD.unit.make_bcc 
# \brief This module is used to generate a unit cell of cscl bcc crystall
# 
# This will work for any crystall type given three primitive vectors of the
# bravais lattice and an origin point. But is designed specifically for a bcc
# crystall

import sys
import os
import numpy as np
import math
from math import sin as sin
from math import cos as cos

#\internal
#
#\brief Generate Rotation Matrix
#\returns rotation matrix for by dergree thetea in 3D
#
#\param theta[thetax,thetay,thetaz]
def rotation(theta):
    tx,ty,tz = theta
    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx),
    cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0,
    cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz),
    0], [0,0,1]])
    return np.dot(Rx, np.dot(Ry, Rz))
## \internal
# \brief Check boundary conditions and remove point if outside
# 
# this is used to get rid of points that fall outside of the box
# I use this because it is simpler to generate too many points and get rid of some
# than it is to generate the right amount
# 
# /returns False if the point falls outside of the simulation box
# /returns x otherwise
# 
# /param x coordinate to check
# /param L box length
def boundary(x,L):
    if x > (L/2.0):
        return False
    if x <= - (L/2.0):
        return False
    if x == 0:
        return True
    return x

## \brief Make a unit crystal from primitive vectors of a cscl crystal
#
# \returns the crystall containing points [atom][type,x,y,z] where type is "A"
# or "B" of the cscl crystall
# \returns an .xyz file called qlattice with the points in vmd type format 
#
# \param a1,a2,a3 are the primitive vectors
# \param basis is a point that will act as the origin of the crystal
# \param L the box size the crystall will fit inside
def make_lattice(a1,a2,a3,basis,L,S=10,name='qlattice.xyz'):
    #MAKE A CRYSTAL
    crystal=[]
    rcrystal=[]
    for i in range(-S,S):
        for j in range(-S,S):
            for k in range(-S,S):
                x = a1[0]*i+a2[0]*j+a3[0]*k+basis[0]
                y = a1[1]*i+a2[1]*j+a3[1]*k+basis[1]
                z = a1[2]*i+a2[2]*j+a3[2]*k+basis[2]
                #apply boundary conditions
                boundary(y,L[1])
                boundary(z,L[2])
                if (boundary(x,L[0])  == False or
                    boundary(y,L[1]) == False or
                    boundary(z,L[2]) == False):
                    pass
                else:
                    #determine if A or B type by checking
                    #for even or odd i +j + k
                    if (i+j+k)%2 ==0:
                        crystal.append([x,y,z,'W'])
                        rcrystal.append([x,y,z])
                    else:
                        crystal.append([x,y,z,'W'])
                        rcrystal.append([x,y,z])
    #write out xyz file for vmd
    fid = open(name,'w')
    fid.write(('%i\n')%(len(crystal)))
    fid.write(('Length L=%.4f\n'%L[0]))
    for i in range(len(crystal)):
        if crystal[i][3] == 'W':
                fid.write(('%c   %f   %f   %f\n')%(crystal[i][3], crystal[i][0],
                crystal[i][1], crystal[i][2]))
    for i in range(len(crystal)):
        if crystal[i][3] == 'B':
                fid.write(('%c   %f   %f   %f\n')%(crystal[i][3], crystal[i][0],
                crystal[i][1], crystal[i][2]))
    fid.close()
    return np.array(rcrystal), crystal



if __name__ == '__main__':
    #sc
    #a1 = np.array([1,0,0])
    #a2 = np.array([0,1,0])
    #a3 = np.array([0,0,1])
    #bcc
    a = 17
    a1 = np.array([-0.5*a,0.5*a,0.5*a])
    a2 = np.array([0.5*a,-0.5*a,0.5*a])
    a3 = np.array([0.5*a,0.5*a,-0.5*a])
    #C0 cyrstal
    #a=9
    #q = 1.1
    #s = 0.3679
    #a1 = np.array([2**(1-1./2*q),2**(1-1./2*q),0])*a
    #a2 = np.array([0,0,2])*a
    #a3 = np.array([-2*s,2*(s+2**(-1./2*q)),1])*a
    #print a1,a2,a3
    ##C1 cyrstal
    #a=9
    #q = 1
    #s = 1/(3*2.**0.5)
    ##a1 = np.array([2**(1-1./2*q),2**(1-1./2*q),0])*a
    ##a2 = np.array([2**(1-1./2*q),0,2**(1-1./2*q)])*a
    ##a3 = np.array([2*(s+2**(-1./2*q)),2*s,-2*s])*a
    #basis = np.array([0.0,0.0,0.0])
    #L=2*np.array([5*a,5*a,5*a])
    #30_20_20 sp8
    #b1 = np.array([-53.75,20.309,15.52])
    #b2 = np.array([-34.25,19.74,14.57])
    #p1 = np.array([-42.9,14.1,-1.19])
    #p2 = np.array([-45.9,3.85,14.45])
    #p3 = np.array([-17.84,29.5,10.74])
    #a1 = p1-b1
    #a2 = p2-b1
    #a3 = p3-b1
    ##
    ##
    #b1 = np.array([-31.299,-10.92,-31.5])
    #p1 = np.array([-24.6,-22.94,-20.88])
    #p2 = np.array([-37.18,-1.2,-18.5])
    #p3 = np.array([-20.23,-5.35,-20.43])
    ##
    ##
    #b1 = np.array([0.0,0.0,0.0])
    #p1 = np.array([-10,0.978,-12.89])
    #p2 = np.array([-9.38,-10.402,-4.156])
    #p3 = np.array([-14.47,2.485,0.98])
   # 
    #a1 = p1-b1
    #a2 = p2-b1
    #a3 = p3-b1
    #b1 = np.array([10.21,-2.57,7.7])
    #basis = b1
    ################
    #roate a bcc lattice
    #################
    #a = 6.
    #print a
    #b1 = np.array([0.0,0.0,0.0])
    #a1 = np.array([0.5*a,0.5*a,0.5*a])
    #a2 = np.array([0.5*a,-0.5*a,0.5*a])
    #a3 = np.array([a,0.0,0.0])
    #R = rotation([math.pi/6,math.pi/6,math.pi/4])
    ##R = rotation([math.pi/6,math.pi/6,math.pi/4])
    #a1 = np.dot(R,a1)
    #a2 = np.dot(R,a2)
    #a3 = np.dot(R,a3)
    #basis = b1
    #count = 0
    #L = [4*a,4*a,4*a]
    #b1 = np.array([-47.209,11.95,10.29])
    #a1 = np.array([-12.773,-6.761,-2.848])
    #a2 = np.array([1.012,-6.126,-13.41])
    #a3 = np.array([-10.056,3.799,-12.633])
    #L = [170*3,170*3,170*3]
    #basis = b1

    VW, VW_names = make_bcc(a1,a2,a3,basis,L,S=30,name='bcc_unit_large.xyz')







