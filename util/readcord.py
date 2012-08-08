import sys
import os
import copy
import numpy as np
import MD.util.readdcd as readdcd
###########################################################
#parse the names from .mol2 file as well as num of particles
############################################################
def name_parse(num_particles,topology):
    #grab the number of particles and move past the garbage
    fid = open(topology,'r')
    fid.readline()
    fid.readline()
    fid.readline()
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
    def __init__(self,topology = 'dna.mol2',trajectory='dna.dcd'):
        self.trajectory=trajectory
        #lets get the number of particles
        #as well as names of particles
        self.DCD = readdcd.DCD(trajectory)
        self.frames = self.DCD.numframes
        self.box_length = self.DCD.box_length_x
        self.box_length_y = self.DCD.box_length_y
        self.box_length_z = self.DCD.box_length_z
        self.num_particles = self.DCD.numatoms
        self.universe,self.num_particles=name_parse(self.num_particles,topology)
        #some info for parsing xyz file
    ########################
    #get the box size at each frame
    #######################
    def box_volume(self):
        self.V = self.DCD.getboxsize()
        return self.V
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
    def matrix_xyz(self,cord=[],cols=3,auto=True,num_frames=1,num_cords=1):
        #frames,rows,columns
        if auto==True:
            self.M = np.zeros((self.num_frames,self.num_cord(cord),cols))
        else:
            self.M = np.zeros((num_frames,num_cords,cols))
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
    #################
    # get index of array coordingates
    #################
    def get_names(self,cord):
        '''grab the indexes of a certain cord'''
        index=[]
        for i,j in enumerate(self.universe):
            for c in cord:
                if j==c:
                    index.append(c)
        return index
    #######################################################################
    # search for all of the lines that have a Cord and store the coordinates
    #######################################################################
    def get_cord(self,cord=[],last=1,start=0,delta=1,get_index='True'):
        if get_index:
            print 'getting index'
            index = self.get_index(cord)
            print index
        else:
            print 'getting universe'
            index = range(len(self.universe))
        print 'finding positions'
        for i,frame in enumerate(range(start,last,delta)):
            f = self.DCD.getframe(frame)
            count=0
            for j in index:
                self.M[i][count][0]=f[j][0]
                self.M[i][count][1]=f[j][1]
                self.M[i][count][2]=f[j][2]
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
        self.num_frames = self.frames
        self.matrix_xyz(cord=cord)
        self.get_cord(cord=cord,last=self.frames)
        return self.return_cords()
    ###########################################################
    #auto cord, note cord must be in the form ['X'] and not 'X'
    #also more than one cord can be passes ie ['X','G'..etc]
    ###########################################################
    def cord_range(self, cord, start=0, delta=1, last=True):
        if last:
            last=self.frames
        self.num_frames = (last-start)/delta
        print 'number of frames'
        print self.num_frames
        print 'Getting Cords',cord
        self.matrix_xyz(cord)
        self.get_cord(cord,last,start=start,delta=delta)
        return self.return_cords()

    ###########################################################
    #auto cord, note cord must be in the form ['X'] and not 'X'
    #also more than one cord can be passes ie ['X','G'..etc]
    ###########################################################
    def all_cord_range(self, start=0, delta=1, last=1):
        self.num_frames = (last-start)/delta
        print 'number of frames'
        print self.num_frames
        print 'Getting Cords'
        self.matrix_xyz(auto=False,num_frames=self.num_frames,num_cords=self.num_particles)
        self.get_cord(last=last,start=start,delta=delta,get_index=False)
        return self.return_cords()

