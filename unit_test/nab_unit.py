import math
import os
import sys
import numpy as np
sys.path.append('/home/cdknorow/Dropbox/Software/MD/base')
import points
def points_on_sphere(N):
    pts = []
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / float(N)
    for k in range(0, N):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return pts

V = points_on_sphere(10000)
phi = 0.575
theta = 3.15
index = []
r =1.0

exclude = 0
for i in range(len(V)):
    if V[i][2] > r*math.cos(phi):
        exclude = 1
    if V[i][2] < -r*math.cos(phi):
        exclude = 1
    if V[i][0] > r*math.cos(phi):
        exclude = 1
    if V[i][0] < -r*math.cos(phi):
        exclude = 1
    if V[i][1] > r*math.cos(phi):
        exclude = 1
    if V[i][1] < -r*math.cos(phi):
        exclude = 1
    if exclude == 1:
        pass
    else:
        index.append(i)
    exclude = 0
Z = []
for i in index:
    Z.append(V[i])
Z = np.array(Z)
#a = 1./3
#b = 2./3
#Z = np.array([[b,b,-a],
#              [a,b,-b],
#              [-a,b,-b],
#              [-b,b,-a],
#              [b,b,a],
#              [a,b,b],
#              [-a,b,b],
#              [-b,b,a],
#              [b,a,-b],
#              [-b,a,-b],
#              [-b,a,b],
#              [b,a,b]])
#Z = np.array([[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]])
#for i in range(len(V)):
#    exclude = 1
#    for j in range(Z.shape[0]):
#        if points.dist(V[i],Z[j],[20,20,20])[0] < 0.1:
#            print 'include'
#            exclude = 0
#    if exclude:
#        pass
#    else:
#        print 'add'
#        index.append(i)
#print index
Z = []
for i in index:
    Z.append(V[i])
Z = np.array(Z)
print Z
#a = 1./3
#b = 2./3
#Z = np.array([[b,b,-a],
#              [a,b,-b],
#              [-a,b,-b],
#              [-b,b,-a],
#              [b,b,a],
#              [a,b,b],
#              [-a,b,b],
#              [-b,b,a]])
fid = open('sphere_points.xyz','w')
fid.write('%i\n\n'%(len(Z)+1))
fid.write('Z 0 0 0\n')
for i in range(Z.shape[0]):
    fid.write('N %.2f %.2f %.2f\n'%(Z[i][0],Z[i][1],Z[i][2]))
fid.close()
u_vectors =[]
#north south and east poles
#for i in range(5):
#    u_vectors.append(np.array(([0, math.sin(i*math.pi/(6*4)),
#        math.cos(i*math.pi/(6.*4))])))
u_vectors.append(np.array(([-1, 0, 0])))
#u_vectors.append(np.array(([0, 1, 0])))
#u_vectors.append(np.array(([0, 0, 1])))
#u_vectors.append(np.array(([1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
##u_vectors.append(np.array(([-1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
###u_vectors.append(np.array(([0, 1/2.**0.5, 1/2.**0.5])))
#u_vectors.append(np.array(([1/2.**0.5, 0, 1/2.**0.5])))
#u_vectors.append(np.array(([1/2.**0.5, 1/2.**0.5, 0])))
#u_vectors.append(np.array(([-1/2.**0.5, 1/2.**0.5, 0])))
#u_vectors.append(np.array(([-1/2.**0.5, 0, 1/2.**0.5])))

####################################################
#Create Master Vector
U = [[] for i in range(len(u_vectors))]
for i in range(len(u_vectors)):
    for j in range(Z.shape[0]):
        if  math.acos(np.dot(u_vectors[i],Z[j])) < theta:
                U[i].append(Z[j])

print len(U[0])
fid = open('u1.xyz','w')
fid.write('%i\n\n'%len(U[0]))
for i in range(len(U[0])):
    fid.write('N %.2f %.2f %.2f\n'%(U[0][i][0],U[0][i][1],U[0][i][2]))
fid.close()
####################################################
#Create Master Vector

def order_parameter(e):
    #bakos 
    b_lamda1 = (13+3*17**0.5)/36.
    b_lamda2 = 5/18.
    b_lamda3 = (13-3*17**0.5)/36.
    #isotropic
    i_lamda1 = 0.701
    i_lamda2 = 0.15
    i_lamda3 = 0.15
    e.sort()
    e = e[::-1]
    print e
    i_order = ((abs(e[0] - i_lamda1) +
             abs(e[1] - i_lamda2) +
             abs(e[2] - i_lamda3)) )
             #(abs(1-i_lamda1)+i_lamda2 + i_lamda3))
    b_order = ((abs(e[0] - b_lamda1) +
             abs(e[1] - b_lamda2) +
             abs(e[2] - b_lamda3)) )
            #(abs(1-b_lamda1)+b_lamda2 + b_lamda3))
    #o = (abs(e[1]-i_lamda2)/(b_lamda2-i_lamda2)  +
    #    abs(e[2]-i_lamda3)/(i_lamda3-b_lamda3)  )
    return i_order/0.6
def Nab(u):
    Quu = np.zeros((3,3))
    for N in range(len(u)):
        for i in range(3):
            for j in range(3):
                Quu[i,j] +=  u[N][i]*u[N][j]
    Quu = Quu / (len(u))
    e_u = np.linalg.eig(Quu)
    return Quu,e_u
for vv,u in enumerate(U):
    Nuu, e_u = Nab(u)
    print 'vectors'
    print u_vectors[vv] 
    print 'Nab'
    print Nuu
    print 'eigenvalues'
    print e_u[0]
