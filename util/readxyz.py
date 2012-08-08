import sys
import os
import copy
import numpy as np
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
###########################################################
#parse the names from .mol2 file as well as num of particles
############################################################
def name_parse_no_mol(trajectory):
    #grab the number of particles and move past the garbage
    fid = open(trajectory,'r')
    num_particles = int(fid.readline().split()[0])
    fid.readline()
    particles=[]
    #get the Letter Name of each particle
    for i in range(num_particles):
        particles.append(fid.readline().split()[0])
    fid.close()
    return particles,num_particles
class ReadCord():
    #data file for all data to get
    def __init__(self,topology='xyz.xyz.xyz',trajectory='xyz.txt',frames=1):
        self.frames = frames
        self.trajectory=trajectory
        #lets get the number of particles
        #as well as names of particles
        try:
            self.universe,self.num_particles=name_parse(topology)
        except:
            self.universe,self.num_particles=name_parse_no_mol(trajectory)
        #some info for parsing xyz file
        fid = open(trajectory,'r')
        self.vmd_step=len(fid.readline())+len(fid.readline())
        self.xyz_step=len(fid.readline())
        fid.close()
        self.frame_step=self.xyz_step*self.num_particles+self.vmd_step
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
    ########################################################################
    #get coorinates of every point in simulation
    ######################################################################
    def get_allcord(self):
        '''this guy isn't implemented anymore never really used it'''
        pass
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
    def get_cord(self,cord):
        index = self.get_index(cord)
        fid = open(self.trajectory,'r')
        for frame in range(self.frames):
            fid.seek(self.vmd_step)
            count=0
            for j in index:
                fid.seek(self.vmd_step+j*self.xyz_step+frame*self.frame_step)
                row = fid.readline().split()
                self.M[frame][count][0]=float(row[1])
                self.M[frame][count][1]=float(row[2])
                self.M[frame][count][2]=float(row[3])
                count+=1
    #######################################################################
    # search for all of the lines that have a Cord and store the coordinates
    #######################################################################
    def get_cord_basic(self,cord):
        index = self.get_index(cord)
        fid = open(self.trajectory,'r')
        for frame in range(self.frames):
            fid.readline()
            fid.readline()
            count=0
            for j in range(self.num_particles):
                row = fid.readline().split()
                if row[0] in cord:
                    self.M[frame][count][0]=float(row[1])
                    self.M[frame][count][1]=float(row[2])
                    self.M[frame][count][2]=float(row[3])
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
    def cord_auto(self,cord):
        print 'Getting Cords',cord
        self.matrix_xyz(cord)
        self.get_cord_basic(cord)
        return self.return_cords()
