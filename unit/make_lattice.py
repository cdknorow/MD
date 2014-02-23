#generate an xml script for an initial configuration
#use coordinates, bonds, type and name data
#use xyz coordinates for the center of the particles
import pickle
import numpy as np
import copy
#zero the dummy particle coordinates
def zero(data):
    num_cords = data[1][0].shape[0]
    center = copy.deepcopy(data[1][0][0])
    print center
    for i in range(data[1].shape[1]):
        data[1][0][i] = data[1][0][i] - center
        print 'this was center'
    for i in range(len(data[3])):
        data[3][i][1] = data[3][i][1] - 28*num_cords
        data[3][i][2] = data[3][i][2] - 28*num_cords
    fid2 = open('test.xyz','w')
    fid2.write('%i\n\n'%data[1].shape[1])
    count = 0
    for i in data[1][0]:
        fid2.write(('%s %.2f %.2f %.2f\n'%(data[0][count],i[0],i[1],i[2])))
        count+=1
    fid2.close()
    return data
#remove positions that are not wanted 
#def slice(VW):
#    count = 0
#    for i in VW:
#        if i[1] <= 37.:
#            if i[2] > -37:
#                if i[0] < 37:
#                    count+=1
#    new_VW = np.zeros((count,3))
#    count = 0
#    for i in VW:
#        if i[1] <37.:
#            if i[2] > -37:
#                if i[0] < 37:
#                    new_VW[count] = i
#                    count += 1
#    return new_VW
#remove positions that are not wanted 
def slice(VW):
    count = 0
    for i in VW:
        if i[0] < 24.5:
                    count+=1
    new_VW = np.zeros((count,3))
    count = 0
    for i in VW:
        if i[0] < 24.5:
            new_VW[count] = i
            count += 1
    return new_VW
#slide the particles in the box and implement images
class image():
    def __init__(self,L):
        self.image = []
        self.L = np.array([L,L,L])
    def slide(self,p):
        self.image.append([0,0,0])
        L = self.L
        if p[0] > self.L[0]/2:
            p[0] = p[0] -L[0]
            self.image[-1][0] = 1
        elif p[0] <- self.L[0]/2:
            p[0] = p[0] + L[0]
            self.image[-1][0] = -1
        if p[1] > self.L[1]/2:
            p[1] = p[1] -L[1]
            self.image[-1][1] = 1
        elif p[1] < -self.L[1]/2:
            p[1] = p[1] + L[1]
            self.image[-1][1] = -1
        if p[2] > self.L[2]/2:
            p[2] = p[2] -L[2]
            self.image[-1][2] = 1
        elif p[2] < -self.L[2]/2:
            p[2] = p[2] + L[2]
            self.image[-1][2] = -1
        return p[0],p[1],p[2]
#write the posision of the particles
def write_position(fid,points, cords,Ntype,L):
    fid2 = open('xyz.xyz','w')
    fid.write('<position>\n')
    fid2.write('%i\n\n'%(len(cords)*len(points)))
    max_z = 0
    max_y = 0
    max_x = 0
    I = image(L)
    count = 0
    for p in points:
        for i in cords:
            x = p[0]+i[0]
            y = p[1]+i[1]
            z = p[2]+i[2]
            x,y,z = I.slide([x,y,z])
            fid.write(('%.2f %.2f %.2f\n'%(x,y,z)))
            fid2.write(('%s %.2f %.2f %.2f\n'%(Ntype[count],x,y,z)))
            count+=1
            if count == len(cords):
                count = 0
            if abs(x) > max_x:
                max_x = abs(x)
            if abs(y) > max_y:
                max_y = abs(y)
            if abs(z) > max_z:
                max_z = abs(z)
    fid2.close()
    print max_x, max_y, max_z
    fid.write('</position>\n')
    fid.write('<image>\n')
    for i in I.image:
        fid.write(('%i %i %i\n'%(i[0],i[1],i[2])))
    fid.write('</image>\n')
#write the body type of particles
def write_body(fid,body,N):
    fid.write('<body>\n')
    for k in range(N):
        for i in body:
            if i == -1:
                fid.write('%i\n'%(i))
            else:
                fid.write('%i\n'%(k))
    fid.write('</body>\n')
#write the Name type of particles
def write_type(fid,Ntype,N):
    fid.write('<type>\n')
    for k in range(N):
        for i in Ntype:
            fid.write('%s\n'%(i))
    fid.write('</type>\n')
#write the bonds
def write_bonds(fid,Bonds,N,num_particle):
    fid.write('<bond>\n')
    for k in range(N):
        for i in Bonds:
            fid.write('%s %i %i\n'%(i[0],i[1]+k*num_particle,i[2]+k*num_particle))
    fid.write('</bond>\n')
if __name__ == "__main__":
    VW = pickle.load(open('bcc_unit.pkl','r'))
    data = pickle.load(open('splice.pkl','r'))
    data = zero(data)
    VW  = slice(VW)
    #write xml nonsense
    N = VW.shape[0]
    num_cords = data[1][0].shape[0]
    L=50.0
    print N
    print num_cords
    fid = open('start.xml','w')
    fid.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    fid.write('<hoomd_xml version="1.4">')
    fid.write('<configuration time_step="0" dimensions="3" natoms=%i">\n'%(N*num_cords))
    fid.write('<box lx="%.2f" ly="%.2f" lz=%.2f" />\n'%(L,L,L))
    write_position(fid,VW,data[1][0],data[0],L)
    write_body(fid,data[2][0],N)
    write_type(fid,data[0],N)
    write_bonds(fid,data[3],N,num_cords)
    fid.write('</configuration>')
    fid.write('</hoomd_xml>')
    fid.close()




