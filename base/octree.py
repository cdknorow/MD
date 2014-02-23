import os
import sys
import numpy as np
import MD.base.points as points


class node():
    def __init__(self,L,delta=4):
        #initialize the box
        delta = L[0]/round(delta)
        self.box = []
        self.L = L
        for i in np.arange(-L[0]/2.,L[0]/2.,delta):
            for j in np.arange(-L[1]/2.,L[1]/2.,delta):
                for k in np.arange(-L[2]/2.,L[2]/2.,delta):
                    self.box.append([np.array([i,j,k]),[],[]])
        #add the box neighbors
        for i in range(len(self.box)):
            for j in range(len(self.box)):
                if points.dist(self.box[i][0],self.box[j][0],L)[0] <= 3**0.5*delta:
                    self.box[i][1].append(j)
    #reset teh box particles
    def clear_particles(self):
        for i in range(self.box):
            self.box[i][2] = []
    #add the particles to their boxes
    def add_particles(self,M):
        self.particle_list = np.zeros((M.shape[0]))
        for i in range(M.shape[0]):
            dsmall = 100
            index = 0
            for j in range(len(self.box)):
                d = points.dist(M[i],self.box[j][0],self.L)[0]
                if d < dsmall:
                    index = j
                    dsmall = d
            self.box[index][2].append(i)
            self.particle_list[i] = int(index)
    #return the index of the particles in neighbor boxes
    def get_particle_neighbors(self,index):
        box_index = int(self.particle_list[index])
        p_index = []
        for i in self.box[box_index][1]:
            p_index.extend(self.box[i][2])
        return p_index
if __name__ == '__main__':
    L = np.array([100,100,100])
    N = node(L,5)



