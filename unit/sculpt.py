import os
import sys

fid = open('thing.xyz','w')


R = 3
total = 27
fid.write('%i \n\n'%total)
fid.write('V 0 0 0 \n')
#edge
E = R +2.5
#corner
C = R+2.0
#center
Z = R+3.5

center =  [[0,0,Z], [0,0,-Z], [0,Z,0], [0,-Z,0],[Z,0,0],[-Z,0,0]]
edge = [[E,0,E],[E,E,0],[0,E,E],[-E,0,-E],[-E,-E,0],[0,-E,-E],
        [E,-E,0],[-E,E,0],[E,0,-E],[-E,0,E],[0,-E,E],[0,E,-E]]
corner = [[C,C,C],[-C,C,C],[C,-C,C],[C,C,-C],
            [-C,-C,C],[C,-C,-C],[-C,C,-C],[-C,-C,-C]]

for i in corner:
    fid.write('Z %.2f %.2f %.2f\n'%(i[0],i[1],i[2]))
for i in center:
    fid.write('Z %.2f %.2f %.2f\n'%(i[0],i[1],i[2]))
for i in edge:
    fid.write('Z %.2f %.2f %.2f\n'%(i[0],i[1],i[2]))


