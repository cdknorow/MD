import os
import sys

fid = open('mylog','w')
fid.write('temperature pressure potential_energy volume\n')
for i in range(20):
    for j in range(1000):
        fid.write('%.2f %.2f %.2f %.2f\n'%(1.0, i, i, i))

