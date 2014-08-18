import os
import sys
import math
sys.path.append('/home/cdknorow/Dropbox/Software/')
import MD


M=MD.ReadCord()
L = M.box_volume()
L_cont = L
delta = 1
fid = open('rotations_split.txt','r')
M = fid.readlines()
fid.close()
x = []
R = []
for line in M:
    s = line.split()
    x.append(int(s[0]))
    R.append(float(s[1]))

out = open('rdiffusion.dat','w')
#sp4
#D = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.06,0.02,0.005,0.002,0.0008,0.00001,
#        0.0001,0.0001,0.0001,0.0001,0.00001,0.00001]
#sp12
#D =     [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#        0.06,0.059,0.04,0.02,0.015,0.01,0.008,0.004,0.003,
#        .0025,0.003,0.003,0.003,0.004,0.004,0.004,0.004,
#        0.004,0.004,0.004]
#nc
#D = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.09,0.1,.01,0.00175,0.0002,0.0002,0.0001,
#        0.00001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.00001,0.00001]

#sp4r41
#D =     [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#        0.07,0.06,0.05,0.04,0.025,0.02,0.0095,0.007,0.004,
#        .0009,0.0005,0.0001]
#sp12
#D  =  [1,1,0.3
#        0.08,0.05,0.03,0.03,0.0325,0.022,0.0245,0.025,0.025,0.015,0.016,
#        0.01,0.011,0.008,0.006,0.004,0.003,0.0025,0.002,0.0015,0.001,
#        0.0009,0.0005,0.0001]
#sp4_n250_sim17
#D  =  [.8,.7,0.6,0.3,0.2,0.08,0.03]
#sp8_sim4
#D = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.9,
#    0.8,0.2,0.2,0.025,0.015,0.008,0.001,0.007,0.007,0.002,0.001,0.005,0.005,0.005,0.004,0.003,0.003,0.002,0.002,0.002,0.002,0.002,0.002,0.001,0.0009,0.0009,0.0009,0.0009,0.0004,0.0004,0.0001,0.0001,0.0001,0.0001]
D =  [1.0,1.0,1.0,1.0,1.0,1.0,1.0,
        0.1, 0.03,0.006,0.001,0.0001,0.001,0.0001,
        0.0001,0.0001,0.01,0.1,0.03,
        0.09,0.009,0.1,0.008,0.003,
        0.001,0.0005,0.0001,0.0001]
count = 0
print len(x)
print len(D)
x = [0]
count = 1
lcut = 5
L_length = [[L_cont[0][0]]]
L_last = L_cont[0]
for i in range(1,len(L_cont),delta):
    if L_cont[i] != L_last :
        L_length[-1].append(count)
        L_length.append([L_cont[i][0]])
        count = 1
        L_last = L_cont[i]
    elif i >= len(L_cont)-delta:
        count+=1
        L_length[-1].append(count)
        break
    else:
        count+=1
count = 0
print L_length
del L_length[-1]
for i in L_length:
    if i[1] > lcut:
        x.append(i[1]+count)
        count+=i[1]
    else:
        print i[1]
        x[-1] += i[1]
        count+=i[1]
print x
t = 0
count = 0
fid2 = open('rcoef.dat','w')
fid2.write('#time L diffusion\n')
for i in range(1,len(x)):
    print i
    for k in range(x[i-1],x[i]):
        if count < len(D):
            out.write(('%i %.2f\n')%(k,math.exp(-2*D[count]*t)))
            t+=1
    if count < len(D):
        fid2.write('%i %.2f %f\n'%(x[i-1],L[x[i-1]][0],D[count]))
        t = 0
        L_last=L[x[i-1]]
        count += 1
        print count,k
out.close()
fid2.close()
import shutil
shutil.copy('rdiffusion.dat','/home/cdknorow/Dropbox')
#t = 0
#count = 0
#fid2 = open('rcoef.dat','w')
#fid2.write('#time L diffusion\n')
#for k in range(len(x)):
#    if count < len(D):
#        out.write(('%i %.2f\n')%(x[k],math.exp(-2*D[count]*t)))
#        if L[k]!=L_last:
#            fid2.write('%i %.2f %f\n'%(x[k],L[k][0],D[count]))
#            t = 0
#            L_last=L[k]
#            count += 1
#            print count,k
#    t+=1
#out.close()
#fid2.close()
#import shutil
#shutil.copy('rdiffusion.dat','/home/cdknorow/Dropbox')
