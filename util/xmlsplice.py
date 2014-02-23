

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
            print i
            print linecache.getline(trajectory,xml[i][0]).split('<')
            s = linecache.getline(trajectory,xml[i][0]).split('<')[1].split()[0]
            xml[i][1] = int(linecache.getline(trajectory,xml[i][0]).split('"')[1])+xml[i][0]+1
            xml.append([xml[i][1]+1,1])
            xmlP[s] = xml[i]
        self.xmlParse = xmlP
    #######################################################################
    # Create x ,y and z plus vectors of each of the proper length as well as distance storage
    #######################################################################
    def matrix_xyz(self,cord='A',cols=3,default = True, args = []):
        #frames,rows,columns
        if default == True:
            self.M = np.zeros((self.frames,self.num_cord(cord),cols))
        else:
            self.M = np.zeros((args['frames'],args['num_cord'],args['cols']))
    #######################################
    # supply the index you want to grab
    ##########################################
    def get_listed_index(self, index, parse):
        count = 0
        for j in index:
            row = linecache.getline(self.trajectory,self.xmlParse[parse][0]+j+1).split()
            for i in range(self.M.shape[2]):
                self.M[0][count][i]=float(row[i])
            count+=1
    def get_listed_type(self, index):
        self.M = []
        for j in index:
            self.M.append(self.universe[j])
    def get_bonds(self,index,parse):
        self.M = []
        for j in range(self.xmlParse[parse][1]-self.xmlParse[parse][0]-1):
            row = linecache.getline(self.trajectory,self.xmlParse[parse][0]+j+1).split()
            if int(row[1]) in index:
                if len(row) == 3:
                    self.M.append([row[0],int(row[1]),int(row[2])])
                if len(row) == 4:
                    self.M.append([row[0],int(row[1]),int(row[2]),int(row[3])])
    ###################################################################
    # return the coordinates
    #################################################################
    def return_cords(self):
        return copy.deepcopy(self.M)
    ###########################################################
    #auto cord, note cord must be in the form ['X'] and not 'X'
    #also more than one cord can be passes ie ['X','G'..etc]
    ###########################################################
    def cord_auto(self,index,parse,arg):
        self.matrix_xyz(default=False,args=arg)
        print self.M.shape
        self.get_listed_index(index,parse)
        return self.return_cords()
    def bonds_auto(self,index,parse,arg):
        self.get_bonds(index,parse)
        return self.return_cords()
    def type_auto(self,index):
        self.get_listed_type(index)
        return self.return_cords()
#splice the particle out
if __name__=='__main__':
    R = ReadCord(trajectory='atoms.dump.0120799935.xml')
    #Index To splice
    index = range(715,1430)
    index = range(20020,20735)
    arg = {'frames':1,'num_cord':len(index),'cols':3}
    P = R.cord_auto(index,'position',arg)
    arg = {'frames':1,'num_cord':len(index),'cols':1}
    B = R.cord_auto(index,'body',arg)
    arg = {'frames':1,'num_cord':len(index),'cols':1}
    T = R.type_auto(index)
    arg = {'frames':1,'num_cord':len(index),'cols':1}
    BO = R.bonds_auto(index,'bond',arg)
    import pickle
    pickle.dump([T,P,B,BO],open('splice.pkl','w'))

