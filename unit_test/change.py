import os 
import sys
import random

fid = open('anticubatic.xyz','r')
out = open('newantic.xyz','w')

M = fid.readlines()
for line in M:
    if 'N' in line:
        if random.random() < 0.5:
            out.write(line.replace('N','S'))
    else:
        out.write(line)
