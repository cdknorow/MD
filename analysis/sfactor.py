# \package MD.analysis.sfactor
# \brief This module finds the static structure factor peaks for a group of
# atoms averaged over several frames.
#
# \b Example
# \code
# stmp,qx = sf.sfactor(VW[start:finish],L=40,l=10) 
# S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
# \endcode
#
# S is the S(q), Q is mag(q_vec), and primitive vectors are q vectors

import os
import numpy as np
import math
import copy

## \brief finds the structure factor.
#
# \returns stmp, qx 
# 
# \param B Matrix of points [frames][atoms][x,y,z]
# \param L periodic boundaries of box 
# \param l number of q_vectors to look for
#
def sfactor(B,L,l=10):
    def do_row(q,r):
        def dot(y,x):
            s=x[:,0]*y[0]+x[:,1]*y[1]+x[:,2]*y[2]
            return s
        a = np.zeros((q.shape[0],1),dtype='complex')
        for i in range(q.shape[0]):
            a[i]=sum(np.exp(1j*dot(q[i],r)))
        a/=r.shape[0]
        return a[:,0]
    #q vectors
    M=np.arange(-l,l+1)
    qx = 2*math.pi/L[0]*M
    qy = 2*math.pi/L[1]*M
    qz = 2*math.pi/L[2]*M
    #q vector matricies
    l=2*l+1
    q=np.zeros((l,3))
    stmp=np.zeros((l,l,l),dtype='complex')
    #Find S(q) for all values of q
    for f in range(B.shape[0]):
        for i in range(l):
            for j in range(l):
                q[:,0] = qx[i]
                q[:,1] = qy[j]
                q[:,2] = qz
                #Sum over all points for given q
                stmp[i,j,:]+=do_row(q,B[f][:])/B.shape[0]
    return stmp,qx


## \brief Sort the 3D peaks of the strucure factor
#
# \returns S(q), mag(q), q vectors
# 
# \param stmp 
# \param qx
#
def sort_sfpeaks(stmp,qx):
    S=[]
    Q=[]
    l = len(qx)
    primitive_vectors=[]
    for i in range(l):
        for j in range(l):
            for k in range(l):
                #magnitude of q vector
                q=(qx[i]*qx[i]+qx[j]*qx[j]+qx[k]*qx[k])**0.5
                #Take the conjugate of stmp because I haven't dont that yet
                S.append((stmp[i][j][k].conjugate()*stmp[i][j][k]).real)
                primitive_vectors.append([qx[i],qx[j],qx[k]])
                #keep track of S and q components
                Q.append(q)
    return S,Q,primitive_vectors


## \brief Filters out S(q) values which are below a certain level.
#
# \return SS is sf peaks with the noise removed
# \return QQ is the corresponding q-vectors magnitude that go with the sf peaks
# \returnPM is the correspoinding q-vector
#
# \param Q |q|
# \param S s(q)
# \praam primitive_vectors q vectors
def sf_filter(Q,S,primitive_vectors,filters=0.05,save='sf'):
    SS=[]
    QQ=[]
    PM=[]
    fid = open(save + '.txt','w')
    for i in range(len(S)):
        if S[i] > filters:
            SS.append(S[i])
            QQ.append(Q[i])
            PM.append(primitive_vectors[i])
    WQ = copy.deepcopy(QQ)
    WQ.sort()
    for i in range(len(WQ)):
        fid.write('%f\n'%WQ[i])
    fid.close()
    return QQ,SS,PM

## \brief Find where sf values should be for bcc.
#
# \returns bars x value, height of bars.
# 
# \param QQ sorted |q| vectors
# \param qpeaks number of bcc peaks to find
def sf_max(QQ,qpeaks=8):
    QQ.sort()
    qstar = QQ[1]
    print 'q* chosen to be'
    print qstar
    height=[]
    bars=[]
    for i in range(1,qpeaks):
        q = qstar*i**0.5
        bars.append(q)
        height.append(1.0)
    return bars,height

## \brief Filters out S(q) values which are below a certain level.
#
# \return SS is sf peaks with the noise removed
# \return QQ is the corresponding q-vectors magnitude that go with the sf peaks
# \returnPM is the correspoinding q-vector
#
# \param Q |q|
# \param S s(q)
# \praam primitive_vectors q vectors
def sf_filter_max(S,Q,primitive_vectors,filters=0.05,save='sf'):
    SS=[]
    QQ=[]
    PM=[]
    for i in range(len(S)):
        if S[i] < filters:
            SS.append(S[i])
            QQ.append(Q[i])
            PM.append(primitive_vectors[i])
    return SS,QQ,PM
## \brief finds the reciprical lattice as well as basis vectors
# 
# \returns realspace, basis vectors
# \returns writes data to file as well
#
# \param QQ input basis vector for q
# \param PM primitive vectors 
#
def recipricol(QQ,PM):
    #delete the 0 value from Q
    try:
        d=QQ.index(0)
        del(QQ[d])
        del(PM[d])
    except:
        print 'no zero value in q'
    #Find the indexes of the min values in QQ that correspond to the first peaks
    #Then add the PM index to b
    b=[]
    q=[]
    fid = open('qnumber.txt','w')
    #record qnumber degeneracy in order
    print 'QQ'
    print QQ
    PP=copy.deepcopy(QQ)
    PP.sort()
    for i in PP:
        fid.write(('%f \n')%(i))
    print QQ
    print 'sorting QQ'
    #remove degeneracy in PM
    for i in range(12):
        ZZ=copy.deepcopy(QQ)
        ZZ.sort()
        index = QQ.index(min(ZZ))
        print min(ZZ)
        b.append(PM[index])
        q.append(QQ[index])
        del(QQ[index])
        del(PM[index])
    print 'finished sorting QQ'
    print len(b)
    #save the basis vectors
    fid = open('latticevec.xyz','w')
    fid.write('13\n')
    fid.write('Atoms\n')
    fid.write('A 0 0 0\n')
    for i in range(12):
        fid.write(('%c   %f   %f   %f\n')%('B',b[i][0],b[i][1], b[i][2]))
    fid.close()
    #find three basis vectors that aren't the same
    print np.array(b)
    b1 = b[int(raw_input('enter value 1?  '))]
    b2 = b[int(raw_input('enter value 2?  '))]
    b3 = b[int(raw_input('enter value 3?  '))]

    print "Found b1, b2, b3!"
    print b1,"\n", b2,"\n",b3
    den=np.dot(b1,np.cross(b2,b3))
    #Find the real-space vectors
    a1 = np.cross(b2,b3)
    a2 = np.cross(b3,b1)
    a3 = np.cross(b1,b2)
    a1 /= den
    a2 /= den
    a3 /= den
    a1 /= 1/(math.pi*2)
    a2 /= 1/(math.pi*2)
    a3 /= 1/(math.pi*2)
    #write to file
    fid =open('recipricol.txt','w')
    fid.write(('b vetcors\n%f %f %f \n %f %f %f \n %f %f %f\n')%(b1[0],
        b1[1],b1[2],b2[0],b2[1],b2[2],b3[0],b3[1],b3[2]))
    fid.write(('a vectors\n%f %f %f \n %f %f %f \n %f %f %f\n')%(a1[0],
        a1[1],a1[2],a2[0],a2[1],a2[2],a3[0],a3[1],a3[2]))
    fid.close()
    b=[b1,b2,b3]
    a=[a1,a2,a3]

    return a, b
## \brief Calculate the debye-waller 
#
# \returns x,y values for debye-waller factor
#
# \params msdr input msdr value
def dbfactor(msdr=5):
    x = np.array(np.arange(0,1.6,0.1))
    y = x*x
    y/=-3
    y*=msdr
    for i in range(len(y)):
        y[i]=math.exp(y[i])
    return x, y
