
import sys
import os
import copy
import numpy as np
import linecache
###########################################################
#parse the names from .mol2 file as well as num of particles
############################################################
def name_parse(topology):
    #grab the number of particles and move past the garbage
    fid = open(topology,'r')
    fid.readline()
    fid.readline()
    num_particles = int(fid.readline().split()[0])
    #more grabage
    fid.readline()
    fid.readline()
    particles=[]
    #get the Letter Name of each particle
    for i in range(num_particles):
        particles.append(fid.readline().split()[1])
    fid.close()
    return particles,num_particles
class ReadCord():
    #data file for all data to get
    def __init__(self,topology='dna.mol2',trajectory='atoms.dump.0003300000.xml',frames=1):
        self.frames = frames
        self.trajectory=trajectory
        #lets get the number of particles
        #as well as names of particles
        try:
            self.universe,self.num_particles=name_parse(topology)
        except:
            self.universe,self.num_particles=name_parse_no_mol(trajectory)
        #generate the info to parse the xml file
        xmlP = {'position':[5,0]}
        xml = [[5,0]]
        total = 5
        for i in range(10):
            s = linecache.getline(trajectory,xml[i][0]).split('<')[1].split()[0]
            xml[i][1] = int(linecache.getline(trajectory,xml[i][0]).split('"')[1])+xml[i][0]+1
            xml.append([xml[i][1]+1,1])
            xmlP[s] = xml[i]
        self.xmlParse = xmlP
    #######################################################################
    # Count all of the lines that have a 'Cord' and return number of counts 
    #######################################################################
    def num_cord(self,cord):
        count=0
        for c in cord:
            count+=self.universe.count(c)
        return count
    #######################################################################
    # Create x ,y and z plus vectors of each of the proper length as well as distance storage
    #######################################################################
    def matrix_xyz(self,cord,cols=3):
        #frames,rows,columns
        self.M = np.zeros((self.frames,self.num_cord(cord),cols))
    #################
    # get index of array coordingates
    #################
    def get_index(self,cord):
        '''grab the indexes of a certain cord'''
        index=[]
        for i,j in enumerate(self.universe):
            for c in cord:
                if j==c:
                    index.append(i)
        index.sort()
        return index
    #######################################################################
    # search for all of the lines that have a Cord and store the coordinates
    #######################################################################
    def get_cord(self,cord,parse):
        index = self.get_index(cord)
        count=0
        for j in index:
            row = linecache.getline(self.trajectory,self.xmlParse[parse][0]+j).split()
            self.M[0][count][0]=float(row[0])
            self.M[0][count][1]=float(row[1])
            self.M[0][count][2]=float(row[2])
            count+=1
    ###################################################################
    # return the coordinates
    #################################################################
    def return_cords(self):
        return copy.deepcopy(self.M)
    ###########################################################
    #auto cord, note cord must be in the form ['X'] and not 'X'
    #also more than one cord can be passes ie ['X','G'..etc]
    ###########################################################
    def cord_auto(self,cord,parse):
        print 'Getting Cords',cord
        self.matrix_xyz(cord)
        print self.M.shape
        self.get_cord(cord,parse)
        return self.return_cords()
if __name__=='__main__':
    R = ReadCord(trajectory='atoms.dump.0075799959.xml')
    M = R.cord_auto('N','acceleration')
    print M
    vx = 0
    vy = 0
    vz = 0
    for i in range(M.shape[1]):
        vx+= M[0][i][0]*M[0][i][0]
        vy+= M[0][i][1]*M[0][i][1]
        vz+= M[0][i][2]*M[0][i][2]
    print vx/M.shape[1],vy/M.shape[1],vz/M.shape[1]
