## \package MD.analysis.cubatic_order
# \brief This scirpt implements an alogrith to find uniaxial and biaxial ordering
#of cubes as defined in ref: http://dx.doi.org/10.1063/1.1711594
import sys
import numpy as np
import MD.base.points as points
# M [cube][x,y,z]
# Z [edge][x,y,z], edge 1-6 represent first cube
def cube_vectors(V,E,L):
    vec = []
    for i in E:
        vec.append(points.unit(points.vector1d(V,i,L)))
    return vec
#kronecker delta function
def delta(i,j):
    if i == j:
        return 1
    else:
        return 0
# \brief find uniaxial and biaxial ordering of cubes within a crystal
#
# \return Q00 and Q22 bond ordering parameter as defined in ref 
# http://dx.doi.org/10.1063/1.1711594
#
# \param M[frame][atoms][x,y,z]: center atom of cube
# \param Z[frame][side][x,y,z]: center points for sides of cube 
# \param L length of box
# \param select reference cube for definiing intial x,y and z axis
def cubatic_order(M,Z,L,select=0):
    #Choose references cube
    ref = M[0][select]
    us = points.unit(points.vector1d(ref,Z[0][select*6],L))
    ws = points.unit(points.vector1d(ref,Z[0][select*6+1],L))
    vs = points.unit(points.vector1d(ref,Z[0][select*6+4],L))
    ####################################################
    #Create u,w,v based on best alignment with reference cube
    u = []
    v = []
    w = []
    for k in range(M.shape[0]):
        for i in range(M.shape[1]):
            if i != select:
                vc = cube_vectors(M[k][i],Z[k][i*6:i*6+6],L)
                usmall = 0
                vsmall = 0
                wsmall = 0
                for j in range(len(vc)):
                    if points.angle_between(us,vc[j])<=points.angle_between(us,vc[usmall]):
                        usmall = j
                    if points.angle_between(ws,vc[j])<=points.angle_between(ws,vc[wsmall]):
                        wsmall = j
                    if points.angle_between(vs,vc[j])<=points.angle_between(vs,vc[vsmall]):
                        vsmall = j
                u.append(vc[usmall])
                v.append(vc[vsmall])
                w.append(vc[wsmall])
            else:
                u.append(ws)
                v.append(us)
                w.append(vs)
    #Find Ordering Tesnors
    Quu = np.zeros((3,3))
    Qww = np.zeros((3,3))
    Qvv = np.zeros((3,3))
    for N in range(len(u)):
        for i in range(3):
            for j in range(3):
                Quu[i,j] += 3. * u[N][i]*u[N][j] - delta(i,j)
                Qvv[i,j] += 3. * v[N][i]*v[N][j] - delta(i,j)
                Qww[i,j] += 3. * w[N][i]*w[N][j] - delta(i,j)
    Quu = Quu / (2.*len(u))
    Qvv = Qvv / (2.*len(v))
    Qww = Qww / (2.*len(w))
    #Find eigenvalues
    e_w = np.linalg.eig(Qww)
    e_u = np.linalg.eig(Quu)
    e_v = np.linalg.eig(Qvv)
    #Identify Laboratory Z,Y and X axis
    def find_max(e):
        evalue = e[0][0]
        index = 0
        if e[0][1] > evalue:
            evalue = e[0][1]
            index = 1
        if e[0][2] > evalue:
            evalue = e[0][2]
            index = 2
        return evalue, index
    #get index of eigenvector and evalue
    e_plus_v = find_max(e_v)
    e_plus_u = find_max(e_u)
    e_plus_w = find_max(e_w)

    # find Z
    s = []
    if e_plus_v[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_v[1][e_plus_v[1]]
        Qz = Qvv
        s.append(0)
    if e_plus_v[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_v[1][e_plus_v[1]]
        Qx = Qvv
        s.append(0)
    if e_plus_u[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_u[1][e_plus_u[1]]
        Qz = Quu
        s.append(1)
    if e_plus_u[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_u[1][e_plus_u[1]]
        Qx = Quu
        s.append(1)
    if e_plus_w[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_w[1][e_plus_w[1]]
        Qz = Qww
        s.append(2)
    if e_plus_w[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_w[1][e_plus_w[1]]
        Qx = Qww
        s.append(2)
    if len(s) != 2:
        print 'error too Qx and Qz were not selected properly'
        asdfasdf

    if 0 not in s:
        y = e_v[1][e_plus_v[1]]
        Qy = Qvv
        s.append(0)
    if 1 not in s:
        y = e_u[1][e_plus_u[1]]
        Qy = Quu
        s.append(1)
    if 2 not in s:
        y = e_w[1][e_plus_w[1]]
        Qy = Qww
        s.append(2)

    x = np.cross(z,y)
    print 'angles between Laboratory'
    print points.angle_between(y,z)

    #Find Q22 and Q00
    Q00 = (np.dot(np.dot(z,Qz),z)+np.dot(np.dot(y,Qy),y)+np.dot(np.dot(x,Qx),x))/3.
    Q22 = (np.dot(np.dot(x,Qx),x) + np.dot(np.dot(y,Qy),y) -
            np.dot(np.dot(x,Qy),x) - np.dot(np.dot(y,Qx),y))/3.

    return Q00, Q22
