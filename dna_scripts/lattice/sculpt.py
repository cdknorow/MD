import os
import sys
import numpy as np
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import MD
import MD.base.points as points

fid = open('thing.xyz','w')


R = 3
total = 27
fid.write('%i \n\n'%total)
fid.write('V 0 0 0 \n')
#edge
E = R +1.95
#corner
C = R+1.35
#center
Z = R+3.0

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

fid.close()

P =  [[0,0,Z], [0,0,-Z], [0,Z,0], [0,-Z,0],[Z,0,0],[-Z,0,0],
        [E,0,E],[E,E,0],[0,E,E],[-E,0,-E],[-E,-E,0],[0,-E,-E],
        [E,-E,0],[-E,E,0],[E,0,-E],[-E,0,E],[0,-E,E],[0,E,-E],
        [C,C,C],[-C,C,C],[C,-C,C],[C,C,-C],[-C,-C,C],[C,-C,-C],
        [-C,C,-C],[-C,-C,-C]]

#get coordinates for Z axis
fid = open('gaussmap_cluster.xyz','r')
fid.readline()
fid.readline()
orientations = []
for line in fid.readlines():
    if line.split()[0] == 'R':
        x = float(line.split()[1])
        y = float(line.split()[2])
        z = float(line.split()[3])
        orientations.append([x,y,z])
fid.close()

#get the coordintes for the unit cell 
fid = open('gaussmap_no_neighbors12.xyz','r')
positions = []
fid.readline()
fid.readline()
for line in fid.readlines():
    if line.split()[0] == 'R':
        x = float(line.split()[1])
        y = float(line.split()[2])
        z = float(line.split()[3])
        positions.append([x,y,z])
fid.close()


#Transform 
L = np.array([50,50,50])
w = np.array([[1,0,0],[0,1,0],[0,0,1]])
V = np.array([0,0,0])
x_r = points.unit(np.array(orientations[1]))
y_r = points.unit(np.array(orientations[2]))
z_r = points.unit(np.array(orientations[3]))
v = np.array([x_r,y_r,z_r])
print points
R = points.reference_rotation(w,v)
location = []
for i in range(len(P)):
    d = points.dist(V,P[i],L)[0]
    c_r = points.unit(points.vector1d(V,P[i],L))
    location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)

#write out 
positions.extend([0,0,0])
fid.open('cube_unit_cell.xyz','w')
fid.write('%i \n\n'%(len(positions)*len(location)))
for i in range(len(positions)):
    for j in range(len(location)):
        x = location[j][0] + positions[i][0]
        y = location[j][1] + positions[i][1]
        z = location[j][2] + positions[i][2]
        fid.write('Z %.2f %.2f %.2f\n'%(x,y,z))
