import os
import sys
import numpy as np
sys.path.append('/home/cdknorow/Dropbox/Software')
import MD
import MD.unit.make_bcc as make_bcc




b1 = np.array([0.0,0.0,0.0])
p1 = np.array([-5.164,12.86,5.1])
p2 = np.array([9.311,10.361,4.205])
p3 = np.array([-4.126,10.339,-9.335])

a1 = p1-b1
a2 = p2-b1
a3 = p3-b1
b1 = np.array([-6.58,-3.91,7.6])
basis = b1
count = 0
for basis in [b1]:
    L = [67,67,67]
    VW, VW_names = make_bcc.make_bcc(a1,a2,a3,basis,L,S=5,name='bcc_unit%i.xyz'%count)
    count+=1
    import pickle
    output = open('bcc_unit.pkl','wb')
    S=pickle.dump(VW,output)
    output.close()
