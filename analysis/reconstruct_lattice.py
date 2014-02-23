import sys
import os
import copy
import numpy as np
import math
## \package MD.analysis.reconstruct_lattice 
# \brief This module is used to reconstruct the real space lattice from a
# structure factor q vectors(note this isn't a very reliable method)
# 
# 

##\ brief Find the real space basis
# 
# This method is primitive and unreliable and will most likley only work for
# simple lattice structures which have only 3 basis vectors
#
# This function also writes out qnumbers.txt which will be an organized list of
# the |q| as well recipricol.txt which holds the recipricol lattice vectors 
#
# \returns real space lattice vectors
# \returns basis vectors in k space
#
# \param QQ |q| where q is the primitive vectors in k space
# \param PM primitive vectors q in k space
def reconstruct_lattice(Q,p,save='recipricol.txt'):
    QQ = copy.deepcopy(Q)
    PM = copy.deepcopy(p)
    #delete the 0 value from Q
    try:
        d=QQ.index(0)
        del(QQ[d])
        del(PM[d])
    except:
        pass
    #Find the indexes of the min values in QQ that correspond to the first
    #peaks of the fcc crystal
    #Then add the PM index to b
    b=[]
    q=[]
    fid = open('qnumber.txt','w')
    PP=copy.deepcopy(QQ)
    PP.sort()
    for i in PP:
        fid.write(('%f \n')%(i))
    for i in range(12):
        ZZ=copy.deepcopy(QQ)
        ZZ.sort()
        index = QQ.index(min(ZZ))
        b.append(PM[index])
        q.append(QQ[index])
        del(QQ[index])
        del(PM[index])
    #Find the R-space vectors
    for i in b:
        print i,"\n"
    sel_1 =int(raw_input('primitive vector 1\n'))
    sel_2 =int(raw_input('primitive vector 2\n'))
    sel_3 =int(raw_input('primitive vector 3\n'))
    b1=b[sel_1]
    b2=b[sel_2]
    b3=b[sel_3]
    print "Found b1,b2,b3"
    fid = open('recip.xyz','w')
    fid.write('37\n\n')
    fid.write('Z 0 0 0\n')
    for i in b[:36]:
        fid.write(('N %f %f %f\n')%(i[0],i[1],i[2]))
    fid.close()
    print b1,"\n", b2,"\n",b3
    fid =open(save,'w')
    fid.write(('%f %f %f %i\n %f %f %f %i\n %f %f %f %i\n')%(b1[0],
        b1[1],b1[2],sel_1,b2[0],b2[1],b2[2],sel_2,b3[0],b3[1],b3[2],sel_3))
    den=np.dot(b1,np.cross(b2,b3))
    a1 = np.cross(b2,b3)
    a2 = np.cross(b3,b1)
    a3 = np.cross(b1,b2)
    a1 /= den
    a2 /= den
    a3 /= den
    a1 /= 1/(math.pi*2)
    a2 /= 1/(math.pi*2)
    a3 /= 1/(math.pi*2)
    fid.write(('%f %f %f \n %f %f %f \n %f %f %f\n')%(a1[0],
        a1[1],a1[2],a2[0],a2[1],a2[2],a3[0],a3[1],a3[2]))
    return a1,a2,a3

