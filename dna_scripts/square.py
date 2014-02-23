import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import util.readxyz as readxyz
import MD.base.points as points
import MD.plot.pyplot as pyplot
import MD.base.clusters as clusters
from MD.analysis.particle_distance import particle_distance 
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import count_neighbors_index
from MD.analysis.nearest_neighbor import second_nearest_neighbors_index
from MD.analysis.nearest_neighbor import min_distance
from MD.plot.histogram import histogram
from MD.plot.histogram import histogram_normal
from MD.plot.histogram import histogram_reg
from MD.analysis.rotation import single_rotation
from MD.analysis.rotation import rotation_angle
from MD.analysis.rotation import diffusion
############ 
# Run from inside a folder
#########
reload(MD)
reload(pyplot)
#########################################
# What you want to do in those directories
########################################
############################################
# returns the average position over a few frames
def Average(VW,L,n_start,n_frames,write=True,save='star.xyz'):
    #Average_position
    Average = np.zeros(VW[0].shape)
    for k in range(n_start,n_start+n_frames):
        for i in range(VW.shape[1]):
            #check the distance between
            for j in range(3):
                if abs(VW[k][i][j]-VW[n_start][i][j]) < L[k][j]/2:
                    Average[i][j] += VW[k][i][j]
                elif VW[k][i][j]-VW[n_start][i][j] < -L[k][j]/2:
                    Average[i][j] += VW[k][i][j]+L[k][j]
                elif VW[k][i][j]-VW[n_start][i][j] > L[k][j]/2:
                    Average[i][j] += VW[k][i][j]-L[k][j]
   # fix any points that may be outside of the box after averaging
    for i in range(Average.shape[0]):
        for j in range(3):
            Average[i][j] /= n_frames
            if Average[i][j] > L[n_start][j]:
                Average[i][j] -= L[n_start][j]
            if Average[i][j] < -L[n_start][j]:
                Average[i][j] += L[n_start][j]
    if write:
        fid = open(save,'w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()
    return Average
#######
# if a point has moved to far we stop using it to find the average position
def Average_simple(VW,L,n_start,n_frames,write=True,cut=10,save='star.xyz'):
    #Average_position
    Average = np.zeros(VW[0].shape)
    for k in range(n_start,n_start+n_frames):
        for i in range(VW.shape[1]):
            #check the distance between
            for j in range(3):
                if abs(VW[k][i][j]-VW[n_start][i][j]) < cut:
                    print abs(VW[k][i][j]-VW[n_start][i][j])
                    Average[i][j] += VW[k][i][j]
   # fix any points that may be outside of the box after averaging
    for i in range(Average.shape[0]):
        for j in range(3):
            Average[i][j] /= n_frames
            if Average[i][j] > L[n_start][j]:
                Average[i][j] -= L[n_start][j]
            if Average[i][j] < -L[n_start][j]:
                Average[i][j] += L[n_start][j]
    if write:
        fid = open(save,'w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()
    return Average
########################
#find the average of particles that are a distance rcut from a point
def Radial_Average(VW,P,L,n_start,n_frames,write=True,save='radial_star.xyz',r_cut=30):
    radial = []
    for i in range(VW.shape[1]):
        d = points.dist(VW[n_start][i],P,L[n_start])
        if d[0] < r_cut:
            radial.append(i)
    VW_R = np.zeros((n_frames,len(radial),3))
    L_R = np.zeros((n_frames,3))
    for k in range(n_frames):
        for i,j in enumerate(radial):
            VW_R[k][i] = VW[k+n_start][j]
            L_R[k] = L[k+n_start]
    return Average_simple(VW_R,L_R,0,n_frames,save=save)
############################################
## msd
############################################
def msd(VW,L,step=1):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x,msd=msd(VW,L,step=step)
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save='MSDtime')
    util.pickle_dump(msd,'msd.pkl')
#get connections
def connections(M):
    try:
        con = util.pickle_load('conn.pkl')
    except:
        fid = open('conn_who.dat','r')
        c_index = M.get_index(['C'])
        g_index = M.get_index(['G'])
        con = []
        count = 0
        for line in fid.readlines():
            print count
            con_K = []
            con_T = []
            #linker alternates between c g and cc gg
            line = line.replace('}','')
            for i in line.split('{')[1].split():
                #translate connections
                try:
                    con_K.append(c_index.index(int(i)))
                except:
                    print 'error connection was not in list'
            for i in line.split('{')[2].split():
                #translate connections
                try:
                    con_T.append(g_index.index(int(i)))
                except:
                    print 'error connection was not in list'
            con.append([con_K,con_T])
            count += 1
        del con[-1]
        util.pickle_dump(con,'conn.pkl')
        return con
############################################
## \brief Find average <u(o)*u(t)> = exp(2Dt)
#
#  \returns diffusion coefficent base on phi
#
# \param Z center of cube sides
# \param VW center of particle
# \param L length of box
def rotation_diffusion_split(Z,VW,L,t0):
    ####
    ####
    sum_t = [0.0 for i in range(VW.shape[0])]
    for k in range(len(t0)-1):
        for i in range(VW.shape[1]):
            theta = diffusion(VW,Z,L,index=i,t0=t0[k],tfinal=t0[k+1])
            for i in range(t0[k],t0[k+1]):
                sum_t[i] += theta[i-t0[k]]
    for i in range(len(sum_t)):
        sum_t[i] = float(sum_t[i])/(VW.shape[1])
    x=range(len(sum_t))
    out = open('rotations_split.txt','w')
    for i in range(len(sum_t)):
        out.write(('%i %.2f\n'%(x[i],sum_t[i])))
    out.close()
    #write ot the ddiffusion coeffiecnent for each t0
    pyplot.plot(x,sum_t,label='theta',
          save='cube_total_rotation_split%i'%i,showleg=True)
## \brief find the largest cubic a cube belongs too
#
# \returns x,y values of histogram (corners, number at corners))
#
# \param M 
# \param VW
# \param L
# \param ndna number of dna per cube
# \param frame frame to look at
def cubics(M, VW, L,
        ndna,frames=[2,10,30,50,100,150,200,250,280,350,400,500,550,590]):
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    try:
        Z = util.pickle_load('Z.pkl')
    except:
        Z = M.cord_auto(['Z'])
        Z = drift_remove(Z,L)
        util.pickle_dump(Z,'Z.pkl')
    #class to do the heavy lifting
    class cubic:
        def __init__(self, center, v1, v2, v3):
            self.cube = []
            self.cubic = []
            self.center = center
            #create the base cube
            self.v1 = v1
            self.v2 = v2
            self.v3 = v3
            for i in range(-1,2):
                for j in range(-1,2):
                    for k in range(-1,2):
                        x = i * self.v1[0] + j * self.v2[0] + k * self.v3[0]
                        y = i * self.v1[1] + j * self.v2[1] + k * self.v3[1]
                        z = i * self.v1[2] + j * self.v2[2] + k * self.v3[2]
                        self.cube.append([x, y, z])
                        self.cubic.append(0)
        #place the point in the right location in the cubic
        def place(self,v,index):
            num = 0
            enum = -1
            for i,j in enumerate(self.cube):
                #we'll assume the largest projection is the closest match
                m =  points.magnitude(points.projection(np.array(j),v))
                if num < m:
                    num = m
                    enum = i
            if enum > 0:
                self.cubic[enum] = index
        #find the max number of cubes in a side
        def count(self):
            #this is the legend to the boxes surrounding 0,0
            box1 = [0,1,3,4,9,10,12,13]
            box2 = [1,2,4,5,10,11,13,14]
            box3 = [3,4,6,7,12,13,15,16]
            box4 = [4,5,7,8,13,14,16,17]
            box5 = [9,10,12,13,18,19,21,22]
            box6 = [10,11,13,14,19,20,22,23]
            box7 = [12,13,15,16,21,22,24,25]
            box8 = [13,14,16,17,22,23,25,26]

            box = [box1,box2,box3,box4,box5,box6,box7,box8]
            max_count = 0
            for b in box:
                count = 0
                for i in b:
                    if self.cubic[i] > 0:
                        count += 1
                if count > max_count:
                    max_count = count
            return max_count+1
        #write out the corrdinates
        def write(self):
            fid = open('cube.xyz','w')
            fid.write(('%i\n')%(len(self.cubic)))
            fid.write('Atoms\n')
            for i in self.cube:
                fid.write(('E %.2f %.2f %.2f\n')%(i[0],i[1],i[2]))
            fid.close()
            fid = open('vector.xyz','w')
            fid.write('4\n')
            fid.write('Atoms\n')
            fid.write(('S %.2f %.2f %.2f\n')%(0,0,0))
            fid.write(('E %.2f %.2f %.2f\n')%(self.v1[0],
                self.v1[1], self.v1[2]))
            fid.write(('E %.2f %.2f %.2f\n')%(self.v2[0],
                self.v2[1], self.v2[2]))
            fid.write(('E %.2f %.2f %.2f\n')%(self.v3[0],
                self.v3[1], self.v3[2]))
            fid.close()
    #get the connections
    try:
        connections=util.pickle_load('concube.pkl')
    except:
        C=M.cord_auto(['C'])
        G=M.cord_auto(['G'])
        connections = con.connections(C[frames],G[frames],L)
        util.pickle_dump(connections,'concube.pkl')
    #find a graph of the particles
    import MD.analysis.graph as graph
    import networkx as nx
    frame_count = 0
    frame_boxes = []
    for k in frames:
        networks, num_networks, deg, neighbors, num_nieghbors, gr = graph.grapher([connections[frame_count]],VW.shape[1],ndna)
        frame_count += 1
        boxes = []
        for i in range(VW.shape[1]):
            #get a list of the particles which are connected
            NN = gr.adjacency_list()
            #get the base vector for the index we are looking at
            #print points.magnitude(points.vector(VW[k][i],Z[k][i*6],L))
            v1 = points.unit(points.vector(VW[k][i],Z[k][i*6+1],L))
            v2 = points.unit(points.vector(VW[k][i],Z[k][i*6+2],L))
            v3 = points.unit(points.vector(VW[k][i],Z[k][i*6+4],L))
            c = cubic(VW[k][i], v1, v2, v3)
            for j in NN[i]:
                #get the vector to the nearest neighbors
                vnear = points.vector(VW[k][i],VW[k][j],L)
                distance = points.magnitude(vnear)*3/2.
                c.place(points.unit(vnear),j)
                for n in NN[j]:
                     if n != i:
                         #look at the nieghbors of the neighbor
                         vfar = points.vector(VW[k][i],VW[k][n],L)
                         #make sure the point isn't into the cube further away
                         if points.magnitude(vfar) < distance:
                             c.place(points.unit(vfar),n)
            boxes.append(c.count())
        hist_theta, xtheta = histogram_reg(boxes,9)
        pyplot.plot_bar(xtheta,hist_theta,
                xlabel='corners',ylabel='number per frame',
                save='cubics', width=0.05)
        box = [0,0,0,0,0,0,0,0,0]
        for r in boxes:
            box[r] += 1
        frame_boxes.append(box)
    y = [[],[],[],[],[],[],[],[]]
    x = []
    label = []
    for i in frame_boxes:
        for j in range(len(y)):
            y[j].append(i[j])
    for i in range(len(y)):
        label.append('%i'%i)
        x.append(frames)
    print frame_boxes
    pyplot.plot_multi(frames,y,label,
            xlabel='corners',ylabel='number per frame',
            save='cubics_all')
############################################
## \brief Find average <u(o)*u(t)> = exp(2Dt)
#
#  \returns diffusion coefficent base on phi
#
# \param Z center of cube sides
# \param VW center of particle
# \param L length of box
def rotation_diffusion_split_length_set(Z,VW,L,t0):
    ####
    ####
    sum_t = [0.0 for i in range(VW.shape[0])]
    for k in range(len(t0)-1):
        for i in range(VW.shape[1]):
            theta = diffusion(VW,Z,L,index=i,t0=t0[k],tfinal=t0[k+1])
            for i in range(t0[k],t0[k+1]):
                sum_t[i] += theta[i-t0[k]]
    for i in range(len(sum_t)):
        sum_t[i] = float(sum_t[i])/(VW.shape[1])
    x=range(len(sum_t))
    out = open('rotations_split.txt','w')
    for i in range(len(sum_t)):
        out.write(('%i %.2f\n'%(x[i],sum_t[i])))
    out.close()
    #write ot the ddiffusion coeffiecnent for each t0
    pyplot.plot(x,sum_t,label='theta',
          save='cube_total_rotation_split%i'%i,showleg=True)
## \brief read mylog file and makes a plot of the energy
#
# \returns x coordinates of time step
# \returns y whatever value of the row you are looking for
#
# \param row the row you would like to grab data from 
def mylog(row = 3):
    #read log file
    import re
    fid = open('mylog.log','r')
    text = fid.readlines()
    fid.close()
    #loop through and pull out values we want
    x = []
    y = []
    for line in text[1:]:
        x.append(line.split()[0])
        y.append(line.split()[row])
    #make a simple plot
    label = re.sub("_"," ",text[0].split()[row])
    save = text[0].split()[row]
    pyplot.plot(x, y, xlabel='Time', ylabel=label, save=save)
    return x, y
## \brief center to end distance of hybridized vs nonhybridized polymers
#
# \returns average distance from center of particle to end of polymer 
#           for hybridizations and free polymer
#
# \param M - ReadCord Class  
# \param L 
#
# V ----- A C k
# W ----- F G T
def center_end_connected(M,L):
    #find the length of the polymer
    import MD.analysis.connections as con
    try:
        K=util.pickle_load('K.pkl')
        T=util.pickle_load('T.pkl')
        V=util.pickle_load('V.pkl')
        W=util.pickle_load('W.pkl')
    except:
        V=M.cord_auto(['V'])
        W=M.cord_auto(['W'])
        K=M.cord_auto(['K'])
        T=M.cord_auto(['T'])
        util.pickle_dump(K,'K.pkl')
        util.pickle_dump(T,'T.pkl')
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(W,'W.pkl')
    #get connections
    try:
        con = util.pickle_load('conn.pkl')
    except:
        fid = open('conn_who.dat','r')
        c_index = M.get_index(['C'])
        g_index = M.get_index(['G'])
        con = []
        for line in fid.readlines():
            con_K = []
            con_T = []
            line = line.replace('}','')
            for i in line.split('{')[1].split():
                #translate connections
                con_K.append(c_index.index(int(i)))
            for i in line.split('{')[2].split():
                #translate connections
                con_T.append(g_index.index(int(i)))
            con.append([con_K,con_T])
        del con[-1]
        util.pickle_dump(con,'conn.pkl')
        #except:
        #    print 'conn_who not created yet'
    h_sum = []
    f_sum = []
    if len(con) != K.shape[0]:
        print len(con)
        print K.shape
        print 'connections and K have different lenghts'
        asdf
    for k in range(len(con)):
        print k
        hybrid = []
        free = []
        #get the points that are not connected
        for i in range(K.shape[1]):
            if i in con[k][0]:
                hybrid.append(points.dist(K[k][i],V[k][i/(K.shape[1]/V.shape[1])],L[k])[0]+0.5)
            else:
                free.append(points.dist(K[k][i],V[k][i/(K.shape[1]/V.shape[1])],L[k])[0]+0.5)
        for i in range(T.shape[1]):
            if i in con[k][1]:
                hybrid.append(points.dist(T[k][i],W[k][i/(T.shape[1]/W.shape[1])],L[k])[0]+0.5)
            else:
                free.append(points.dist(T[k][i],W[k][i/(T.shape[1]/W.shape[1])],L[k])[0]+0.5)
        h_sum.append(sum(hybrid)/len(hybrid))
        f_sum.append(sum(free)/len(free))
    f_sum[-1]= f_sum[-2]
    h_sum[-1]= h_sum[-2]
    x = range(V.shape[0])
    print 'number of DNA per NP'
    print K.shape[1]/V.shape[1]
    print T.shape[1]/W.shape[1]
    pyplot.plot2(x,h_sum,x,f_sum,xlabel='timestep',ylabel='sigma',label1='hybridizations',
        label2='free',save='free_hybrid_connections',showleg='true')
## \brief end to end distance of hybridized vs nonhybridized polymers
#
# \returns average distance from start to end of polymer 
#           for hybridizations and free polymer
#
# \param M - ReadCord Class  
# \param L 
#
# V ----- A C k
# W ----- F G T
def end_end_connected(M,L): #find the length of the polymer
    import MD.analysis.connections as con
    try:
        K=util.pickle_load('K.pkl')
        T=util.pickle_load('T.pkl')
        S=util.pickle_load('S.pkl')
    except:
        K=M.cord_auto(['K'])
        T=M.cord_auto(['T'])
        S=M.cord_auto(['M'])
        util.pickle_dump(K,'K.pkl')
        util.pickle_dump(T,'T.pkl')
        util.pickle_dump(S,'S.pkl')
    #get connections
    try:
        con = util.pickle_load('conn.pkl')
    except:
        fid = open('conn_who.dat','r')
        c_index = M.get_index(['C'])
        g_index = M.get_index(['G'])
        con = []
        for line in fid.readlines():
            con_K = []
            con_T = []
            line = line.replace('}','')
            for i in line.split('{')[1].split():
                #translate connections
                con_K.append(c_index.index(int(i)))
            for i in line.split('{')[2].split():
                #translate connections
                con_T.append(g_index.index(int(i)))
            con.append([con_K,con_T])
        del con[-1]
        util.pickle_dump(con,'conn.pkl')
        #except:
        #    print 'conn_who not created yet'
    h_sum = []
    f_sum = []
    for k in range(len(con)):
        print k
        hybrid = []
        free = []
        #get the points that are not connected
        for i in range(K.shape[1]):
            if i in con[k][0]:
                hybrid.append(points.dist(K[k][i],S[k][i],L[k])[0]+0.5)
            else:
                free.append(points.dist(K[k][i],S[k][i],L[k])[0]+0.5)
        for i in range(T.shape[1]):
            if i in con[k][1]:
                hybrid.append(points.dist(T[k][i],S[k][i+K.shape[1]],L[k])[0]+0.5)
            else:
                free.append(points.dist(T[k][i],S[k][i+K.shape[1]],L[k])[0]+0.5)
        print hybrid
        h_sum.append(sum(hybrid)/len(hybrid))
        f_sum.append(sum(free)/len(free))
    f_sum[-1]= f_sum[-2]
    h_sum[-1]= h_sum[-2]
    x = range(S.shape[0])
    pyplot.plot2(x,h_sum,x,f_sum,xlabel='timestep',ylabel='sigma',label1='hybridizations',
        label2='free',save='free_hybrid_end_end_connections',showleg='true')
#\brief find the gauss map of connections surrounding NC
def linker_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    try:
        #V = util.pickle_load('V.pkl')
        C = util.pickle_load('C.pkl')
        CC = util.pickle_load('G.pkl')
    except:
         C = M.cord_auto(['C'])
         CC = M.cord_auto(['G'])
         util.pickle_dump(C,'C.pkl')
         util.pickle_dump(CC,'G.pkl')
    try:
        conn = util.pickle_load('conn.pkl')
    except:
        conn = connections(M)
        util.pickle_dump(conn,'conn.pkl')
    A = VW.shape[1]/2
    ndnaA=C.shape[1]/(VW.shape[1]/2)
    ndnaB=CC.shape[1]/(VW.shape[1]/2)
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    for k in frames:
        print k
        conn_location = []
        try:
            for i in conn[k][0]:
                #we must rotate about a specific cube reference frame
                cube = i/ndnaA
                if cube == 1:
                    V = VW[k][cube]
                    x_r = points.unit(points.vector1d(V,Z[k][cube*6+1],L))
                    y_r = points.unit(points.vector1d(V,Z[k][cube*6+2],L))
                    z_r = points.unit(points.vector1d(V,Z[k][cube*6+5],L))
                    v = np.array([x_r,y_r,z_r])
                    R = points.reference_rotation(v,w)
                    d = points.dist(V,C[k][i],L)[0]
                    c_r = points.unit(points.vector1d(V,C[k][i],L))
                    conn_location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        except:
            print 'something went wrong here'
        #for i in conn[k][1]:
        #    #get the correct cube
        #    cube = i/ndnaB + A
        #    V = VW[k][cube]
        #    x_r = points.unit(points.vector1d(V,Z[k][cube*6],L))
        #    y_r = points.unit(points.vector1d(V,Z[k][cube*6+1],L))
        #    z_r = points.unit(points.vector1d(V,Z[k][cube*6+5],L))
        #    v = np.array([x_r,y_r,z_r])
        #    d = points.dist(V,CC[k][i],L)[0]
        #    R = points.reference_rotation(v,w)
        #    c_r = points.unit(points.vector1d(V,CC[k][i],L))
        #    conn_location.append(points.unit(np.dot(R,c_r))*d)
        gauss_map.append(conn_location)
    #########
    fid = open('gaussmap_connections.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
    fid.close()
    ###########
    fid = open('gaussmap_connections_total.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        max_gauss+=len(k)
    fid.write('%i\n\n'%(max_gauss+4))
    fid.write('E  1 0 0\n')
    fid.write('E  0 1 0\n')
    fid.write('E  0 0 1\n')
    fid.write('V  0 0 0\n')
    for k in range(len(gauss_map)):
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    fid.close()
    ############
    fid = open('gaussmap_connections_step.xyz','w')
    max_gauss = 0
    step = 10
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    max_gauss = max_gauss*step
    for s in range(0,len(gauss_map),step):
        fid.write('%i\n\n'%(max_gauss+4))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        count = 4
        for k in range(s,s+step): 
            for i in gauss_map[k]:
                fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
                count+=1
        for i in range(max_gauss+4 - count):
            fid.write('N  0 0 0\n')
    fid.close()
#\brief find the gauss map of polymers surrounding NC
def polymer_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    try:
        #V = util.pickle_load('V.pkl')
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    A = VW.shape[1]/2
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ndna = P.shape[1]/VW.shape[1]
    print P.shape
    for k in frames:
        location = []
        print k
        try:
            for i in range(VW.shape[1]):
                #we must rotate about a specific cube reference frame
                if i in range(1):
                    V = VW[k][i]
                    x_r = points.unit(points.vector1d(V,Z[k][i*6+1],L[k]))
                    y_r = points.unit(points.vector1d(V,Z[k][i*6+2],L[k]))
                    z_r = points.unit(points.vector1d(V,Z[k][i*6+5],L[k]))
                    v = np.array([x_r,y_r,z_r])
                    R = points.reference_rotation(v,w)
                    for j in range(1,ndna,2):
                        d = points.dist(V,P[k][j+i*ndna],L[k])[0]
                        c_r = points.unit(points.vector1d(V,P[k][j+i*ndna],L[k]))
                        location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        except:
            print 'something went wrong here'
        gauss_map.append(location)
    #########
    #fid = open('gaussmap_polymer.xyz','w')
    #max_gauss = 0
    #for k in gauss_map:
    #    if len(k) > max_gauss:
    #        max_gauss = len(k)
    #for k in range(len(gauss_map)):
    #    fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
    #    fid.write('E  1 0 0\n')
    #    fid.write('E  0 1 0\n')
    #    fid.write('E  0 0 1\n')
    #    fid.write('V  0 0 0\n')
    #    for i in gauss_map[k]:
    #        fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #    for i in range(max_gauss - len(gauss_map[k])):
    #       fid.write('N  0 0 0\n')
    #fid.close()
    ###########
    fid = open('gaussmap_polymers_total.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    max_gauss = max_gauss*len(gauss_map)
    fid.write('%i\n\n'%(max_gauss+4))
    fid.write('E  1 0 0\n')
    fid.write('E  0 1 0\n')
    fid.write('E  0 0 1\n')
    fid.write('V  0 0 0\n')
    count = 4
    for k in range(len(gauss_map)):
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            count+=1
    fid.close()
#\brief find the gauss map of polymers surrounding NC
def polymer_inter_gauss(M,VW, Z,L, frames, rcut=20, step=5e4):
    try:
        #V = util.pickle_load('V.pkl')
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    A = VW.shape[1]/2
    gauss_map = []
    gauss_inter = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ndna = P.shape[1]/VW.shape[1]
    min_D = rcut
    print P.shape
    for k in frames:
        location = []
        location_other = []
        print k
        try:
            for i in range(VW.shape[1]):
                #we must rotate about a specific cube reference frame
                if i in range(1):
                    V = VW[k][i]
                    x_r = points.unit(points.vector1d(V,Z[k][i*6+1],L[k]))
                    y_r = points.unit(points.vector1d(V,Z[k][i*6+2],L[k]))
                    z_r = points.unit(points.vector1d(V,Z[k][i*6+5],L[k]))
                    v = np.array([x_r,y_r,z_r])
                    R = points.reference_rotation(w,v)
                    same_set = []
                    for j in range(1,ndna,2):
                        d = points.dist(V,P[k][j+i*ndna],L[k])[0]
                        c_r = points.unit(points.vector1d(V,P[k][j+i*ndna],L[k]))
                        #location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
                        print points.dist(V,P[k][j+1*ndna],L[k])[1]
                        location.append(points.dist(V,P[k][j+i*ndna],L[k])[1])
                        same_set.append(j+i*ndna)
                    for poly in range(1,P.shape[1],2):
                        if poly not in same_set:
                            d = points.dist(V,P[k][poly],L[k])[0]
                            if d < rcut:
                                if d< min_D:
                                    min_D = d
                                c_r = points.unit(points.vector1d(V,P[k][poly],L[k]))
                                #location_other.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
                                location_other.append(points.dist(V,P[k][poly],L[k])[1])
        except:
            print 'ERROR: Failed to rotate frame!'
        gauss_map.append(location)
        gauss_inter.append(location_other)
    #########
    #fid = open('gaussmap_polymer.xyz','w')
    #max_gauss = 0
    #for k in gauss_map:
    #    if len(k) > max_gauss:
    #        max_gauss = len(k)
    #for k in range(len(gauss_map)):
    #    fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
    #    fid.write('E  1 0 0\n')
    #    fid.write('E  0 1 0\n')
    #    fid.write('E  0 0 1\n')
    #    fid.write('V  0 0 0\n')
    #    for i in gauss_map[k]:
    #        fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #    for i in range(max_gauss - len(gauss_map[k])):
    #       fid.write('N  0 0 0\n')
    #fid.close()
    ###########
    fid = open('gaussmap_polymers_total.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    max_gauss = max_gauss*len(gauss_map)
    fid.write('%i\n\n'%(max_gauss+4))
    fid.write('E  1 0 0\n')
    fid.write('E  0 1 0\n')
    fid.write('E  0 0 1\n')
    fid.write('V  0 0 0\n')
    count = 4
    for k in range(len(gauss_map)):
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            count+=1
    fid.close()
    fid = open('gaussmap_polymers_inter.xyz','w')
    max_gauss = 0
    for k in gauss_inter:
        if len(k) > max_gauss:
            max_gauss = len(k)
    max_gauss = max_gauss*len(gauss_inter)
    fid.write('%i\n\n'%(max_gauss))
    fid.write('V  0 0 0\n')
    count = 1
    R = rcut-min_D
    delta = R/4.0 
    for k in range(len(gauss_inter)):
        for i in gauss_inter[k]:
            d = points.dist(np.array([0,0,0]),i,L[frames[k]])[0]
            if d < min_D+delta:
                fid.write('A  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            elif d < min_D+delta*2:
                fid.write('B  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            elif d < min_D+delta*3:
                fid.write('C  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            else:
                fid.write('D  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            count+=1
    for i in range(max_gauss-count):
            fid.write('T  0 0 0\n')
    fid.close()
#\brief find the gauss map of linkers surrounding NC
def polymer_rho(M,L,frames):
    import MD.analysis.polymer as polymer
    reload(polymer)
    B = polymer.analyze_polymer(M,L,frames)
    B.generate_rotation_map()
    hmulti_x1 = []
    hmulti_y1 = []
    hmulti_x2 = []
    hmulti_y2 = []
    hmulti_x3 = []
    hmulti_y3 = []
    hmulti_y4 = []
    hmulti_x4 = []
    hmulti_y5 = []
    hmulti_x5 = []
    hmulti_y6 = []
    hmulti_x6 = []
    hmulti_y7 = []
    hmulti_x7 = []
    label = []
    label_theory = []
    for count in range(len(B.rotation_map)):
        zed,zed_total,e2e,c2e,c2s,curl,theta,theta_ideal,rho,radial =  B.get_polymer_stats(count)
        average_e2e  = sum(e2e)/len(e2e)
        #hist_zed , xe, max_z = histogram_normal(zed,20)
        #pyplot.plot_bar(xe,hist_zed,save='rho_zed_distance',xmax=10,ymax=0.3)
        #plt.close()
        xf = np.arange(0.2,5,0.5)
        hist_frac_e2e = B.avg_radial(rho,e2e,xf)
        hist_frac_zed = B.avg_radial(rho,zed,xf)
        hist_frac_zed_total = B.avg_radial(rho,zed_total,xf)
        hist_frac_c2e = B.avg_radial(rho,c2e,xf)
        hist_frac_theta = B.avg_radial(rho,theta,xf)
        hist_frac_theta_ideal = B.avg_radial(rho,theta_ideal,xf)
        hist_frac_curl = B.avg_radial(rho,curl,xf)
        hmulti_x1.append(xf)
        hmulti_y1.append(hist_frac_e2e)
        hmulti_x2.append(xf)
        hmulti_y2.append(hist_frac_zed)
        hmulti_x3.append(xf)
        hmulti_y3.append(hist_frac_c2e)
        hmulti_x4.append(xf)
        hmulti_y4.append(hist_frac_theta)
        hmulti_x5.append(xf)
        hmulti_y5.append(hist_frac_theta_ideal)
        hmulti_x6.append(xf)
        hmulti_y6.append(hist_frac_curl)
        hmulti_x7.append(xf)
        hmulti_y7.append(hist_frac_zed_total)
        #thoeretical values
        theory_zed = []
        theory_c2e = []
        theory_theta = []
        for i in range(len(xf)):
            Z = (average_e2e) / ((1 + xf[i]**2/B.edge**2)**(0.5))
            theory_zed.append(Z)
        for i in range(len(xf)):
            theory_c2e.append((average_e2e)+(xf[i]**2+B.edge**2)**0.5)
        for i in range(len(xf)):
            t  = math.atan(B.edge/xf[i])
            theory_theta.append(t)
        hmulti_x2.append(xf)
        hmulti_y2.append(theory_zed)
        hmulti_x3.append(xf)
        hmulti_y3.append(theory_c2e)
        hmulti_x4.append(xf)
        hmulti_y4.append(theory_theta)
        hmulti_x7.append(xf)
        hmulti_y7.append(theory_zed)
        label.append('%.2f'%L[frames[count]][0])
        label_theory.append('actual L=%.2f'%L[frames[count]][0])
        label_theory.append('theory L=%.2f'%L[frames[count]][0])
    pyplot.plot_multi(hmulti_x1,hmulti_y1,label,xlabel='rho',ylabel='sigma',
            title='Polymer end to end distance', save='rho_e2e')
    pyplot.plot_multi(hmulti_x2,hmulti_y2,label_theory,xlabel='rho',ylabel='sigma',
            title='Z projection',save='rho_zed',split=True)
    pyplot.plot_multi(hmulti_x3,hmulti_y3,label_theory,xlabel=r'$\rho(\sigma)$',
            ylabel=r'center to polymer end($\sigma$)',
            title= 'Center of Cube to Polymer end',save='rho_c2e',split=True)
    pyplot.plot_multi(hmulti_x4,hmulti_y4,label_theory,xlabel=r'$\rho(\sigma)$',
            ylabel='theta',
            title= 'theta',save='rho_theta',split=True)
    pyplot.plot_multi(hmulti_x5,hmulti_y5,label,xlabel=r'$\rho(\sigma)$',
            ylabel='deviation from ideal theta',
            title= 'theta',save='rho_theta_ideal')
    pyplot.plot_multi(hmulti_x6,hmulti_y6,label,xlabel=r'$\rho(\sigma)$',
            ylabel='curling', title= 'curling',save='rho_curl')
    pyplot.plot_multi(hmulti_x7,hmulti_y7,label_theory,xlabel=r'$\rho(\sigma)$',
            ylabel='azimuth', title= 'top azimuth',save='rho_top_zed',split=True)
    B.write_variables(hmulti_x5,hmulti_y5,label,'rho_theta_ideal.dat')
    B.write_with_theory(hmulti_x7,hmulti_y7,label_theory,'rho_top_zed.dat')
def polymer_radial(M,L,frames):
    import MD.analysis.polymer as polymer
    reload(polymer)
    B = polymer.analyze_polymer(M,L,frames)
    B.generate_rotation_map()
    hmulti_x1 = []
    hmulti_y1 = []
    hmulti_x2 = []
    hmulti_y2 = []
    hmulti_x3 = []
    hmulti_y3 = []
    hmulti_y4 = []
    hmulti_x4 = []
    hmulti_y5 = []
    hmulti_x5 = []
    hmulti_y6 = []
    hmulti_x6 = []
    hmulti_y7 = []
    hmulti_x7 = []
    label = []
    label_theory = []
    for count in range(len(B.rotation_map)):
        zed,zed_total,e2e,c2e,c2s,curl,theta,theta_ideal,rho,radial = B.get_polymer_stats(count)
        average_e2e  = sum(e2e)/len(e2e)
        #hist_zed , xe, max_z = histogram_normal(zed,20)
        #pyplot.plot_bar(xe,hist_zed,save='radial_zed_distance',xmax=10,ymax=0.3)
        #plt.close()
        xf = np.arange(0.2,12,0.5)
        hist_frac_e2e = B.avg_radial(radial,e2e,xf)
        hist_frac_zed = B.avg_radial(radial,zed,xf)
        hist_frac_zed_total = B.avg_radial(radial,zed_total,xf)
        hist_frac_c2e = B.avg_radial(radial,c2e,xf)
        hist_frac_theta = B.avg_radial(radial,theta,xf)
        hist_frac_theta_ideal = B.avg_radial(radial,theta_ideal,xf)
        hist_frac_curl = B.avg_radial(radial,curl,xf)
        hmulti_x1.append(xf)
        hmulti_y1.append(hist_frac_e2e)
        hmulti_x2.append(xf)
        hmulti_y2.append(hist_frac_zed)
        hmulti_x3.append(xf)
        hmulti_y3.append(hist_frac_c2e)
        hmulti_x4.append(xf)
        hmulti_y4.append(hist_frac_theta)
        hmulti_x5.append(xf)
        hmulti_y5.append(hist_frac_theta_ideal)
        hmulti_x6.append(xf)
        hmulti_y6.append(hist_frac_curl)
        hmulti_x7.append(xf)
        hmulti_y7.append(hist_frac_zed_total)
        #thoeretical values
        theory_zed = []
        theory_c2e = []
        theory_theta = []
        theory_x = []
        xf_theory = np.arange(0.2,5,0.5)
        for i in range(len(xf_theory)):
            Z = average_e2e/(1+xf_theory[i]**2/B.edge**2)**0.5
            R = xf_theory[i]/B.edge*average_e2e/(1+xf_theory[i]**2/B.edge**2)**0.5
            t  = math.atan(B.edge/xf_theory[i])
            theory_zed.append(Z)
            theory_x.append(R+xf_theory[i])
            theory_theta.append(t)
        hmulti_x2.append(theory_x)
        hmulti_y2.append(theory_zed)
        hmulti_x4.append(theory_x)
        hmulti_y4.append(theory_theta)
        hmulti_x7.append(theory_x)
        hmulti_y7.append(theory_zed)
        label.append('%.2f'%L[frames[count]][0])
        label_theory.append('actual L=%.2f'%L[frames[count]][0])
        label_theory.append('theory L=%.2f'%L[frames[count]][0])
    pyplot.plot_multi(hmulti_x1,hmulti_y1,label,xlabel='radial',ylabel='sigma',
            title='Polymer end to end distance', save='radial_e2e')
    pyplot.plot_multi(hmulti_x2,hmulti_y2,label_theory,xlabel='radial',ylabel='sigma',
            title='Z projection',save='radial_zed',split=True)
    pyplot.plot_multi(hmulti_x3,hmulti_y3,label,xlabel=r'radial',
            ylabel=r'center to polymer end($\sigma$)',
            title= 'Center of Cube to Polymer end',save='radial_c2e',split=True)
    pyplot.plot_multi(hmulti_x4,hmulti_y4,label_theory,xlabel='radial',
            ylabel='theta',
            title= 'theta',save='radial_theta',split=True)
    pyplot.plot_multi(hmulti_x5,hmulti_y5,label,xlabel='radial',
            ylabel='deviation from ideal theta',
            title= 'theta',save='radial_theta_ideal')
    pyplot.plot_multi(hmulti_x6,hmulti_y6,label,xlabel='radial',
            ylabel='curling', title= 'curling',save='radial_curl')
    pyplot.plot_multi(hmulti_x7,hmulti_y7,label_theory,xlabel='radial',
            ylabel='zed', title= 'top zed',save='radial_top_zed',split=True)
    B.write_with_theory(hmulti_x7,hmulti_y7,label_theory,'radial_top_zed.dat')
#\brief find the gauss map of linkers surrounding NC
def draw_polymer_azimuth(M,L, frames, rcut=1.0, step=5e4):
    print "getting cors"
    P = M.cord_frames(['M'],frames)
    S = M.cord_frames(['S'],frames)
    VW = M.cord_frames(['V','W'],frames)
    Z = M.cord_frames(['Z'],frames)
    print "finished getting cords"
    A = VW.shape[1]/2
    gauss_map = []
    gauss_map_poly = []
    print VW.shape
    print P.shape
    print S.shape
    print Z.shape
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ndna = P.shape[1]/VW.shape[1]
    ndna_poly = S.shape[1]/VW.shape[1]
    spacers = S.shape[1]/VW.shape[1]/(ndna/2)
    number_spacers = 6
    print P.shape
    L_frame=[]
    #distance from center to edge of cube
    edge = points.dist(VW[0][0],Z[0][1],L[0])[0]+0.84
    #draw edges of the cube
    def write_line(fid,c,p1,p2):
        s1 = "{%.2f %.2f %.2f}"%(p1[0],p1[1],p1[2])
        s2 = "{%.2f %.2f %.2f}"%(p2[0],p2[1],p2[2])
        fid.write(("draw line %s %s width 3\n")%(s1,s2))
    for count,k in enumerate(frames):
        L_frame.append(k)
        location = []
        location_poly = []
        print k
        for i in range(VW.shape[1]):
            #we must rotate about a specific cube reference frame
            if i in range(10):
                V = VW[count][i]
                x_r = points.unit(points.vector1d(V,Z[count][i*6+1],L[k]))
                y_r = points.unit(points.vector1d(V,Z[count][i*6+2],L[k]))
                z_r = points.unit(points.vector1d(V,Z[count][i*6+5],L[k]))
                v = np.array([x_r,y_r,z_r])
                R = points.reference_rotation(v,w)
                for j in range(ndna):
                    d = points.dist(V,P[count][j+i*ndna],L[k])[0]
                    c_r = points.unit(points.vector1d(V,P[count][j+i*ndna],L[k]))
                    location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
                for j in range(spacers-1,ndna_poly,spacers):
                    for jj in range(number_spacers):
                        d = points.dist(V,S[count][j-jj+i*ndna_poly],L[k])[0]
                        c_r = points.unit(points.vector1d(V,S[count][j-jj+i*ndna_poly],L[k]))
                        location_poly.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        gauss_map.append(location)
        gauss_map_poly.append(location_poly)
    #find the azimuthal and the end 2 end disance of polymer
    label=[]
    for count,k in enumerate(gauss_map):
        azi = []
        e2e = []
        radi = []
        curl = [ ]
        azi_total = []
        for p in range(0,len(k),2):
            #find the azimuth 
            s = []
            for i in k[p]:
                s.append(abs(i))
            i = s.index(max(s))
            if i == 0:
                radi.append(points.dist(s,np.array([edge,0,0]),L[L_frame[count]])[0])
            if i == 1:
                radi.append(points.dist(s,np.array([0,edge,0]),L[L_frame[count]])[0])
            if i == 2:
                radi.append(points.dist(s,np.array([0,0,edge]),L[L_frame[count]])[0])
            azi.append(abs(k[p+1][i])-edge)
            e2e.append(points.dist(k[p],k[p+1],L[L_frame[count]])[0])
            curler = []
            for j in range(spacers):
                curler.append(abs(k[p+1][i])-abs(gauss_map_poly[count][3*p/2+number_spacers-(j+1)][i]))
            curl.append(min(curler))
            if curl[-1] < 0:
                azi_total.append(azi[-1]+abs(curl[-1]))
                new_theta=True
                i = curler.index(min(curler))
            else:
                new_theta=False
                azi_total.append(azi[-1])
            try:
                if new_theta == True:
                    tot_e2e = points.dist(k[p],gauss_map_poly[count][3*p/2+number_spacers-(1+i)],L[L_frame[count]])[0]
                    theta.append(math.asin(azi_total[-1]/tot_e2e))
                else:
                    theta.append(math.asin(azi_total[-1]/e2e[-1]))
            except:
                theta.append(math.pi/2)
            theta_ideal.append(math.pi-math.acos((c2s[-1]**2+e2e[-1]**2-c2e[-1]**2)/(2*e2e[-1]*c2s[-1])))
def draw_polymer(M,L, frames,NP=[1,2], rcut=1.0, step=5e4):
    print "getting cors"
    P = M.cord_frames(['M'],frames)
    VW = M.cord_frames(['V','W'],frames)
    print "finished getting cords"
    A = VW.shape[1]/2
    ndna = P.shape[1]/VW.shape[1]
    #distance from center to edge of cube
    #draw edges of the cube
    def write_line(fid,p1,p2):
        s1 = "{%.2f %.2f %.2f}"%(p1[0],p1[1],p1[2])
        s2 = "{%.2f %.2f %.2f}"%(p2[0],p2[1],p2[2])
        fid.write(("draw arrow %s %s\n")%(s1,s2))
    c = ['green','blue','red']
    for count,k in enumerate(frames):
        fid = open('polymer%i.tcl'%k,'w')
        fid.write('proc vmd_draw_arrow {mol start end} {\n' +
                  'set middle [vecadd $start [vecscale 0.9 [vecsub $end'+
                  ' $start]]]\n'+'graphics $mol cylinder $start $middle' +
                  ' radius 0.15\n'+'graphics $mol cone $middle $end radius' +
                  ' 0.25\n}\n')
        #for i in range(VW.shape[1]):
        for cc,i in enumerate(NP):
            color = c[cc]
            fid.write('draw color '+color+'\n')
            for j in range(0,ndna,2):
                d = points.dist(P[count][j+i*ndna],P[count][j+1+i*ndna],L[k])[1]
                write_line(fid,P[count][j+i*ndna],P[count][j+i*ndna]+d)
        fid.close()
def polymer_interpenetrate(M,VW, Z,L, frames, rcut=20.0, step=5e4):
    try:
        #V = util.pickle_load('V.pkl')
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    A = VW.shape[1]/2
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ndna = P.shape[1]/VW.shape[1]
    print P.shape
    L_frame=[]
    for k in frames:
        print k
        L_frame.append(k)
        other = []
        same = []
        for i in range(VW.shape[1]):
            if i in range(10):
                V = VW[k][i]
                same_set = []
                for j in range(ndna):
                    d = points.dist(V,P[k][j+i*ndna],L[k])[0]
                    same_set.append(j+i*ndna)
                    same.append(d)
                for poly in range(1,P.shape[1],2):
                    if poly not in same_set:
                        d = points.dist(V,P[k][poly],L[k])[0]
                        if d < rcut:
                            other.append(d)
        gauss_map.append([same,other])
    count = 0
    for k in gauss_map:
        print count
        same = []
        other = []
        hist_same , xs, max_z = histogram_normal(k[0],50)
        hist_other, xo, max_z = histogram_normal(k[1],50)
        pyplot.plot_bar(xs,hist_same,save=False,xmax=20,ymax=0.3,color='r',width=0.5)
        pyplot.plot_bar(xo,hist_other,save='compare_distance%2f'%L[L_frame[count]][0],xmax=20,ymax=0.1,color='b',width=0.5,xlabel='sigma',ylabel='k')
        plt.close()
        #pyplot.plot_bar(xo,hist_other,save='other_distance%2f'%L[L_frame[count]][0],xmax=20,ymax=0.1,color='b',width=0.5)
        #plt.close()
        count+=1
# FInd the rotational motion of cube
#####################
def rotations(VW,Z,L,chord):
    x=range(0,VW.shape[1],10)
    for i in x:
        theta, phi = single_rotation_box(VW,Z,L,index=i)
        x=np.arange(theta.shape[0])
        pyplot.plot2(x,theta,x,phi,label1='theta',
                label2='phi',save='rotation%i'%i,showleg=True)
    sum_t = [0.0 for i in range(VW.shape[0]-1)]
    sum_p = [0.0 for i in range(VW.shape[0]-1)]
    for k in range(VW.shape[1]):
        theta, phi = single_rotation(VW,Z,L,index=k)
        for i in range(len(theta)-1):
            x = abs(theta[i+1]-theta[i])
            if x>math.pi:
                x = abs(x - 2*math.pi)
            sum_t[i] += x
        for i in range(len(phi)-1):
            y = abs(phi[i+1]-phi[i])
            if y>math.pi/2:
                y = abs(y - math.pi)
            sum_p[i] += y
    for i in range(len(sum_t)-1):
        sum_t[i+1] += sum_t[i]
        sum_p[i+1] += sum_p[i]
    for i in range(len(sum_t)):
        sum_t[i] = sum_t[i]/VW.shape[1]
        sum_p[i] = sum_p[i]/VW.shape[1]
    x=range(len(sum_p))
    pyplot.plot2(x,sum_t,x,sum_p,label1='theta',
            label2='phi',save='cube_total_rotation%i'%i,showleg=True)
    return x, theta, phi
#\brief find the angle the polymer makes with the grafting surface
# with respect to its radial distance from the center of the edge
#
# uses law of cosines with
#           a = center to start of polymer
#           b = polymer end to end
#           c = polymer center to end
#           r = rho distance 
def polymer_angle(M,VW,Z, L, frames, rcut=12, step=5e4,delta=5):
    try:
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    ndna = P.shape[1]/VW.shape[1]
    edge = points.dist(VW[0][0],Z[0][0],L[0])[0]
    multi_x = []
    multi_y = []
    label = []
    for k in frames:
        print k
        theta = []
        r = []
        for dt in range(delta):
            try:
                for i in range(VW.shape[1]):
                    #we must rotate about a specific cube reference frame
                    if i in range(VW.shape[1]):
                        V = VW[k+dt][i]
                        for j in range(1,ndna,2):
                            a = points.dist(V,P[k+dt][j+i*ndna-1],L[k+dt])[0]
                            b = points.dist(P[k+dt][j+i*ndna],P[k+dt][j+i*ndna-1],L[k+dt])[0]
                            c = points.dist(V,P[k+dt][j+i*ndna],L[k+dt])[0]
                            theta.append(np.arccos((a**2 + b**2 - c**2) / (2 * a * b))-math.pi/2.)
                            r.append((a**2-edge**2)**0.5)
            except:
                print 'ERROR: Failed to rotate frame!'
            #plot some distance vs its radial distance
        def avg_radial(radi,y,xf):
            hist_frac_s = [[0.0,0] for i in range(len(xf))]
            for i in range(len(radi)):
                r = 0
                for j in range(len(xf)):
                    if radi[i] > xf[j]:
                        r = j
                hist_frac_s[r][0] += y[i]
                hist_frac_s[r][1] += 1
            hist_frac =[]
            for i in hist_frac_s:
                if i[1] != 0:
                    hist_frac.append(i[0]/i[1])
                else:
                    hist_frac.append(i[0])
            return hist_frac
        xf_theta = np.arange(0,6,0.3)
        hist_theta = avg_radial(r,theta,xf_theta)
        multi_x.append(xf_theta)
        multi_y.append(hist_theta)
    #output to text file
    fid = open('theta_multi.dat','w')
    fid.write('#rho')
    for i in range(len(frames)):
        label.append('L=%.1f'%L[frames[i]][0])
        fid.write('  '+label[-1])
    fid.write('\n')
    for i in range(len(multi_x[0])):
        fid.write("  %.2f"%multi_x[0][i]) 
        for j in range(len(multi_y)):
            fid.write("  %.2f"%multi_y[j][i]) 
        fid.write('\n')
    pyplot.plot_multi(multi_x,multi_y,label,xlabel='rho',ylabel='average theta',
            title='Polymer theta', save='theta_multi')
#\brief find the angle the polymer makes with the grafting surface
# with respect to its radial distance from the center of the edge
#
# uses law of cosines with
#           a = center to start of polymer
#           b = polymer end to end
#           c = polymer center to end
#           r = rho distance 
def polymer_angle_hist(M,VW,Z, L, frames, rcut=12, step=5e4, delta=10):
    try:
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    ndna = P.shape[1]/VW.shape[1]
    edge = points.dist(VW[0][0],Z[0][0],L[0])[0]
    multi_x = []
    multi_y = []
    label = []
    for k in frames:
        theta = []
        r = []
        print k
        for dt in range(delta):
            try:
                for i in range(VW.shape[1]):
                    #we must rotate about a specific cube reference frame
                    if i in range(VW.shape[1]):
                        V = VW[k+dt][i]
                        for j in range(1,ndna,2):
                            a = points.dist(V,P[k+dt][j+i*ndna-1],L[k+dt])[0]
                            b = points.dist(P[k+dt][j+i*ndna],P[k+dt][j+i*ndna-1],L[k+dt])[0]
                            c = points.dist(V,P[k+dt][j+i*ndna],L[k+dt])[0]
                            theta.append(np.arccos((a**2 + b**2 - c**2) / (2 * a * b))-math.pi/2.)
                            r.append((a**2-edge**2)**0.5)
            except:
                print 'ERROR: Failed to rotate frame!'
            #plot some distance vs its radial distance
        def split_radial(radi,yf):
            ay = []
            by = []
            cy = []
            for i in range(len(radi)):
                if radi[i] < 3.1:
                    ay.append(yf[i])
                elif radi[i] < 4.1:
                    by.append(yf[i])
                else:
                    cy.append(yf[i])
            ahist, ax, max_z = histogram_normal(ay,20)
            bhist, bx, max_z = histogram_normal(by,20)
            chist, cx, max_z = histogram_normal(cy,20)
            return [ahist,ax],[bhist,bx],[chist,cx]
            return hist_frac
        a,b,c = split_radial(r,theta)
        multi_x = [a[1],b[1],c[1]]
        multi_y = [a[0],b[0],c[0]]
        label = ['rho<3','rho<4','rho<5']
    pyplot.plot_multi(multi_x,multi_y,label,xlabel='rho',ylabel='hist theta',
            title='Polymer theta', save='radial_hist%.1f'%L[frames[0]][0])
def nn_distance(VW, L, frames, rcut=20.0, step=5e4):
    #plot average end to end distance vs radial distance
    def avg_radial(y,xf):
        hist_frac_s = [[0.0,0] for i in range(len(xf))]
        for i in range(len(y)):
            r = 0
            for j in range(len(xf)):
                if y[i] > xf[j]:
                    r = j
            hist_frac_s[r][0] += y[i]
            hist_frac_s[r][1] += 1
        hist_frac =[]
        for i in hist_frac_s:
            if i[1] != 0:
                hist_frac.append(i[0]/i[1])
            else:
                hist_frac.append(i[0])
        return hist_frac
    L_frame=[]
    count = 0
    hmulti_x = []
    hmulti_y = []
    label = []
    counter = []
    for k in frames:
        print k
        L_frame.append(k)
        distance = []
        for i in range(VW.shape[1]):
            for j in range(VW.shape[1]):
                d = points.dist(VW[k][i],VW[k][j],L[k])[0]
                if d < rcut:
                    distance.append(d)
        x = np.arange(5,rcut,0.2)
        if count >8:
            if count%2 == 0:
                hmulti_y.append(avg_radial(distance,x))
                hmulti_x.append(x)
                label.append('%.2f'%L[L_frame[count]][0])
                counter.append(count)
        count+=1
    x_bar = []
    y_bar = []
    for i in range(len(hmulti_y)):
        first = True
        for j in range(len(hmulti_y[i])):
            if hmulti_y[i][j] != 0 and first == True:
                y_bar.append(hmulti_x[i][j])
                x_bar.append(L[L_frame[counter[i]]][0])
                first = False
    pyplot.plot(x_bar,y_bar,xlabel='box length',ylabel='interparticle',
                save='interpart box')
    pyplot.plot_multi(hmulti_x,hmulti_y,label,xlabel='sigma',ylabel='interparticle',
                save='interpart')
def ssDNA_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    try:
        #V = util.pickle_load('V.pkl')
        P = util.pickle_load('SSDNA.pkl')
    except:
         P = M.cord_auto(['K','T'])
         util.pickle_dump(P,'SSDNA.pkl')
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ndna = P.shape[1]/VW.shape[1]
    print P.shape
    for k in frames:
        location = []
        print k
        try:
            for i in range(VW.shape[1]):
                #we must rotate about a specific cube reference frame
                if i in range(20):
                    V = VW[k][i]
                    x_r = points.unit(points.vector1d(V,Z[k][i*6+1],L[k]))
                    y_r = points.unit(points.vector1d(V,Z[k][i*6+2],L[k]))
                    z_r = points.unit(points.vector1d(V,Z[k][i*6+5],L[k]))
                    dd = points.dist(V,Z[k][i*6+5],L[k])[0]
                    v = np.array([x_r,y_r,z_r])
                    R = points.reference_rotation(v,w)
                    for j in range(ndna):
                        d = points.dist(V,P[k][j+i*ndna],L[k])[0]
                        c_r = points.unit(points.vector1d(V,P[k][j+i*ndna],L[k]))
                        location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        except:
            print 'something went wrong here'
        gauss_map.append(location)
    #########
    fid = open('gaussmap_dna.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
        fid.write('E  %.2f 0 0\n'%dd)
        fid.write('E  0 %.2f 0\n'%dd)
        fid.write('E  0 0 %.2f\n'%dd)
        fid.write('V  0 0 0\n')
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
           fid.write('N  0 0 0\n')
    fid.close()
    ###########
    fid = open('gaussmap_dna_step.xyz','w')
    max_gauss = 0
    step = 10
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    max_gauss = max_gauss*step
    for s in range(0,len(gauss_map)-step,step):
        fid.write('%i\n\n'%(max_gauss+4))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        count = 4
        for k in range(s,s+step): 
            for i in gauss_map[k]:
                fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
                count+=1
        for i in range(max_gauss+4 - count):
            fid.write('N  0 0 0\n')
    fid.close()
#\brief find the gauss map of neighbors surrounding NC
# and allign them with the axis of center cube
def neighbor_gauss(VW, Z,L, frames, step=5e4,rcut=25):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    A = VW.shape[1]/2
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    for k in frames:
        print k
        N_location = []
        for i in range(VW.shape[1]):
            if i in range(150):
                #we must rotate about a specific cube reference frame
                x_r = points.unit(points.vector1d(VW[k][i],Z[k][i*6+1],L[k]))
                y_r = points.unit(points.vector1d(VW[k][i],Z[k][i*6+2],L[k]))
                z_r = points.unit(points.vector1d(VW[k][i],Z[k][i*6+5],L[k]))
                v = np.array([x_r,y_r,z_r])
                R = points.reference_rotation(v,w)
                for N in count_neighbors_index(VW[k],i,L[k],count=8,rcut=rcut)[0]:
                    d = points.dist(VW[k][i],VW[k][N],L[k])[0]
                    c_r = points.unit(points.vector1d(VW[k][i],VW[k][N],L[k]))
                    N_location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        gauss_map.append(N_location)
    #########
    fid = open('gaussmap_neighbors.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
    fid.close()
    ############
    #fid = open('gaussmap_neighbors_total.xyz','w')
    #max_gauss = 0
    #for k in gauss_map:
    #    max_gauss+=len(k)
    #fid.write('%i\n\n'%(max_gauss+4))
    #fid.write('E  1 0 0\n')
    #fid.write('E  0 1 0\n')
    #fid.write('E  0 0 1\n')
    #fid.write('V  0 0 0\n')
    #for k in range(len(gauss_map)):
    #    for i in gauss_map[k]:
    #        fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #fid.close()
    ############
    #fid = open('gaussmap_neighbors_step.xyz','w')
    #max_gauss = 0
    #step = 10
    #for k in gauss_map:
    #    if len(k) > max_gauss:
    #        max_gauss = len(k)
    #max_gauss = max_gauss*step
    #for s in range(0,len(gauss_map)-step,step):
    #    fid.write('%i\n\n'%(max_gauss+4))
    #    fid.write('E  1 0 0\n')
    #    fid.write('E  0 1 0\n')
    #    fid.write('E  0 0 1\n')
    #    fid.write('V  0 0 0\n')
    #    count = 4
    #    for k in range(s,s+step):
    #        for i in gauss_map[k]:
    #            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #            count+=1
    #    for i in range(max_gauss+4 - count):
    #        fid.write('N  0 0 0\n')
    #fid.close()
#\brief find the gauss map of neighbors surrounding NC
# but do not allighn with axis of cube
def neighbor_no_gauss(VW, L, frames, step=5e4,rcut=25,count=8):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    gauss_map = []
    print VW.shape
    for k in frames:
        print k
        N_location = []
        #VW_avg = Average(VW,L,k,5,save='star.xyz')
        VW_avg = Radial_Average(VW,VW[k][1010],L,k,5,save='radial_star.xyz')
        print VW_avg
        for i in range(VW_avg.shape[0]):
            #we must rotate about a specific cube reference frame
            counter = 0
            for N in count_neighbors_index(VW_avg,i,L[k],count=count,rcut=rcut)[0]:
                counter+=1
                d = points.dist(VW_avg[i],VW_avg[N],L[k])[1]
                N_location.append(d)
            if counter>0:
                pass
        gauss_map.append(N_location)
    #########
    fid = open(('gaussmap_no_neighbors%i.xyz'%count),'w')
    max_gauss = 0
    CM = []
    for g in gauss_map:
        print '.'
        CM.append(clusters.cluster(g,rcut=1.5,cluster_cut=9))
        index = []
        #sort the cluster
        for i in CM[0]:
            d_small = 100
            j_index = 0
            for j in range(len(CM[-1])):
                d =  points.dist_np(CM[-1][j],i)[0]
                if d < d_small:
                    d_small =  d
                    j_index = j
            index.append(j_index)
        C_step = []
        print index
        for i in index:
            C_step.append(CM[-1][i])
        CM[-1] = C_step
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4+len(CM[k]),frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        #write out center of mass of clusters
        for i in CM[k]:
            fid.write('R  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
    fid.close()
    ############
    #fid = open('gaussmap_no_neighbors_total.xyz','w')
    #max_gauss = 0
    #for k in gauss_map:
    #    max_gauss+=len(k)
    #fid.write('%i\n\n'%(max_gauss+4))
    #fid.write('E  1 0 0\n')
    #fid.write('E  0 1 0\n')
    #fid.write('E  0 0 1\n')
    #fid.write('V  0 0 0\n')
    #for k in range(len(gauss_map)):
    #    for i in gauss_map[k]:
    #        fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #fid.close()
    #############
    #fid = open('gaussmap_no_neighbors_step.xyz','w')
    #max_gauss = 0
    #step = 10
    #for k in gauss_map:
    #    if len(k) > max_gauss:
    #        max_gauss = len(k)
    #max_gauss = max_gauss*step
    #for s in range(0,len(gauss_map)-step,step):
    #    fid.write('%i\n\n'%(max_gauss+4))
    #    fid.write('E  1 0 0\n')
    #    fid.write('E  0 1 0\n')
    #    fid.write('E  0 0 1\n')
    #    fid.write('V  0 0 0\n')
    #    count = 4
    #    for k in range(s,s+step):
    #        for i in gauss_map[k]:
    #            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #            count+=1
    #    for i in range(max_gauss+4 - count):
    #        fid.write('N  0 0 0\n')
    #fid.close()
#\brief find for a single cube
def neighbor_no_gauss_individual(VW, L, frames, step=5e4,rcut=25,count=8):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    A = VW.shape[1]/2
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    for k in frames:
        print k
        N_location = []
        for i in range(VW.shape[1]):
            if i == 1:
                #we must rotate about a specific cube reference frame
                for N in count_neighbors_index(VW[k],i,L[k],count=count,rcut=rcut)[0]:
                    d = points.dist(VW[k][i],VW[k][N],L[k])[1]
                    N_location.append(d)
        gauss_map.append(N_location)
    #########
    fid = open(('gaussmap_no_neighbors%i.xyz'%count),'w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
    fid.close()
    ############
    #fid = open('gaussmap_no_neighbors_total.xyz','w')
    #max_gauss = 0
    #for k in gauss_map:
    #    max_gauss+=len(k)
    #fid.write('%i\n\n'%(max_gauss+4))
    #fid.write('E  1 0 0\n')
    #fid.write('E  0 1 0\n')
    #fid.write('E  0 0 1\n')
    #fid.write('V  0 0 0\n')
    #for k in range(len(gauss_map)):
    #    for i in gauss_map[k]:
    #        fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #fid.close()
    #############
    #fid = open('gaussmap_no_neighbors_step.xyz','w')
    #max_gauss = 0
    #step = 10
    #for k in gauss_map:
    #    if len(k) > max_gauss:
    #        max_gauss = len(k)
    #max_gauss = max_gauss*step
    #for s in range(0,len(gauss_map)-step,step):
    #    fid.write('%i\n\n'%(max_gauss+4))
    #    fid.write('E  1 0 0\n')
    #    fid.write('E  0 1 0\n')
    #    fid.write('E  0 0 1\n')
    #    fid.write('V  0 0 0\n')
    #    count = 4
    #    for k in range(s,s+step):
    #        for i in gauss_map[k]:
    #            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    #            count+=1
    #    for i in range(max_gauss+4 - count):
    #        fid.write('N  0 0 0\n')
    #fid.close()
def linker_gauss_test(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    try:
        N = util.pickle_load('N.pkl')
    except:
        N = M.cord_auto(['N'])
        util.pickle_dump(N,'N.pkl')
    A = VW.shape[1]/2
    nsphere = N.shape[1]/VW.shape[1]
    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    for k in frames:
        print k
        conn_location = []
        for i in range(N.shape[1]):
            #we must rotate about a specific cube reference frame
            cube = i/nsphere
            if cube == 3:
                V = VW[k][cube]
                x_r = points.unit(points.vector1d(V,Z[k][cube*6+1],L))
                y_r = points.unit(points.vector1d(V,Z[k][cube*6+2],L))
                z_r = points.unit(points.vector1d(V,Z[k][cube*6+5],L))
                v = np.array([x_r,y_r,z_r])
                R = points.reference_rotation(v,w)
                d = points.dist(V,N[k][i],L)[0]
                c_r = points.unit(points.vector1d(V,N[k][i],L))
                print '###'
                print v[1]
                print 'rotate'
                print np.dot(R,np.transpose(v[1]))
                conn_location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        gauss_map.append(conn_location)
    #########
    fid = open('gaussmap_connections.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4,frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
    fid.close()
    ###########
    fid = open('gaussmap_connections_total.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        max_gauss+=len(k)
    fid.write('%i\n\n'%(max_gauss+4))
    fid.write('E  1 0 0\n')
    fid.write('E  0 1 0\n')
    fid.write('E  0 0 1\n')
    fid.write('V  0 0 0\n')
    for k in range(len(gauss_map)):
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
    fid.close()
##################################################### #find number of a specified pair of particles that are within r_cut of each
#also finds the number of connections
#####################################################
def find_networks(M, VW, L, n_finish=1, n_start=0, delta=30, rcut=1.0, ndna=25):
    #The total number of frames we are going to look at
    x=np.arange(n_start,n_finish,delta)
    print len(x)
    #Find the number of connections at specific points x
    import MD.canalysis.connections as con
    import MD.analysis.connections as conn
    import MD.analysis.graph as graph
    try:
        connections = util.pickle_load('con.pkl')
    except:
        try:
            C = util.pickle_load('C.pkl')
            G = util.pickle_load('G.pkl')
        except:
            G=M.cord_auto(['G'])
            C=M.cord_auto(['C'])
            util.pickle_dump(C,'C.pkl')
            util.pickle_dump(G,'G.pkl')
        connections = con.connections(C[x],G[x],L,rcut=rcut)
        util.pickle_dump(connections,'con.pkl')
    #plot total number of connections
    num_connections=conn.num_connections(connections,VW.shape[1])
    con_all = []
    for i in range(len(connections)):
        con_all.append(len(connections[i])/2.)
    pyplot.plot(x,num_connections,xlabel='Time',
            ylabel='hyrbid. density', save='connections')
    pyplot.plot(x,con_all,xlabel='Time',
            ylabel='hyrbid', save='connections_all')
    plt.close()
    #get the info
    networks, num_networks, deg, neighbors, num_n, gr = graph.grapher(connections,VW.shape[1],ndna)
    util.pickle_dump(networks,'net.pkl')
    #plot the number of neighbors at each timesteps
    pyplot_eps.plot(x,num_n,xlabel='Time',
            ylabel='num neighbors', save='neighbors')
    plt.close()
    print 'making plot'
    net = []
    for i in networks:
        net.append(len(i))
    pyplot_eps.plot(x,net,xlabel='t',ylabel='networks',save='net')
    label = ['networks','1','2','3','4','5','6','7','8','9']
    pyplot_eps.plot(x,networks,xlabel='t',ylabel='Networks',save='net')
    pyplot.plot_multi(x,deg,label,xlabel='time',ylabel='number',save='con_net')
    return x, networks, connections
## find the z dinswer question a rectional length of polymer ## away from the surface of the cube
def azimuthal_end(M,VW,L,var):
    KT=M.cord_auto(['K','T'])
    S=M.cord_auto(['M'])
    Z = M.cord_auto(['Z'])
    VW=M.cord_auto(['V','W'])
    z=[]
    ### We need to know which z vector to use
    ### We can find the minimum distance to Z bead from M bead and use that one
    ### to identify the side of the cube that the ssDNA is attatched. 
    def min_distance(Z,S,L):
        num=0
        side=0
        smallest=20
        for i in range(6):
            num = points.dist(Z[i],S,L)[0]
            if smallest > num:
                smallest = num
                side = i
        return side
    #loop through and azimithal distance
    for k in range(VW.shape[0]-5,VW.shape[0]):
        for i in range(var['nsphere']):
            for j in range(var['ndna']):
                #project vector from center to end onto z vector
                side = min_distance(Z[k][i*6:(i+1)*6],S[k][j+i*var['ndna']],L)
                v1=points.vector(VW[k][i],Z[k][side+i*6],L)
                v2=points.vector(VW[k][i],KT[k][j+i*var['ndna']],L)
                proj = points.projection(v1,v2)
                z.append(points.magnitude(proj)-points.magnitude(v1))
    hist_z,xz,max_z=histogram(z,50)
    pyplot.plot_bar(xz,hist_z,save='z_distance')
#################################################
# Finds uniaxial and biaxial ordering 
#################################################
def cubic_order(VW,Z,L_cont,x):
    import MD.analysis.cubatic_order as c
    Q = []
    reload(c)
    frame = []
    print VW.shape
    print len(L_cont)
    xvalue = []
    for k in x:
        #a = x[i+1]-x[i]
        #delta =  a/4
        #k = x[i]+a/2
        #print k
        #print L_cont[k]
        xvalue.append(k)
        #Q.append(c.cubatic_order_eigenvalues(VW[k:k+delta],Z[k:k+delta],
        #        L_cont[k],k,delta=delta))
        delta = 0
        Q.append(c.cubatic_order_eigenvalues(np.array([VW[k]]),np.array([Z[k]]),
                L_cont[k],k,delta=1))
        #if i == 6:
        #    print L_cont[k]
        #    asdfasdfa
        #c.cube_allign(VW[k],Z[k],L_cont[k])
        #c.least_fit_cubes(VW[k],Z[k],L_cont[k],k)
        ##q,q2 = c.cubatic_order(VW[i:i+5],Z[i:i+5],L)
        #Q.append(qsave)
        #Q2.append(q2save)
        #frame.append(i)
    #Q =  Q/max(Q)
    #print Q2
    x = xvalue
    fid = open('cubatic.txt','w')
    for i in range(len(Q)):
        fid.write(('%i %.2f %.2f\n'%(x[i],Q[i],L_cont[x[i]][0])))
    fid.close()

    pyplot.plot(x,Q,save='Qnew')
    #pyplot.plot(x,Q2,save='Q2')
def cubic_order_animate(VW,Z,L_cont,x):
    import MD.analysis.cubatic_order as c
    Q = []
    reload(c)
    frame = []
    x = range(50,100)
    xvalue = []
    for k in x:
        c.least_fit_cubes(VW[k],Z[k],L_cont[k],k)
def find_allignment(VW,Z,L_cont,k):
    import MD.analysis.cubatic_order as c
    Q = []
    reload(c)
    frame = []
    Q.append(c.new_cubatic_order(VW[k],Z[k],
            L_cont[k]))
    fid = open('cubatic.txt','w')
    for i in range(len(Q)):
        fid.write(('%i %.2f\n'%(x[i],Q[i])))
    fid.close()

    pyplot.plot(x,Q,save='Qnew')
    #pyplot.plot(x,Q2,save='Q2')
def print_box_volume(L_cont,delta=10):
    fid = open('box_length.txt','w')
    fid.write('frame  L\n') 
    for i in range(0,len(L_cont)):
        fid.write(('%i %.2f %.2f %.2f\n')%(i,L_cont[i][0],L_cont[i][1],L_cont[i][2]))
    fid.close()
def assymetric_pressure(VW,L,frame):
    fid = open('animate/bakos_index%i.txt'%frame,'r')
    M = []
    for line in  fid.readlines():
        M.append(int(line))
    fid.close()
    import linecache
    pressure1 = linecache.getline('pressure.log',frame+2)
    pp = 0
    pxx= 0
    pyy= 0
    pzz= 0
    for i in range(2,6):
        pressure = linecache.getline('pressure.log',frame+i)
        pp += float(pressure.split()[1])/4.
        pxx += float(pressure.split()[2])/4.
        pyy += float(pressure.split()[3])/4.
        pzz += float(pressure.split()[4])/4.

    xx = []
    yy = []
    zz = []
    #find the 2nd nearest neighbors
    for i in range(VW.shape[1]):
        N =  count_neighbors_index(VW[frame],i,L[frame],count=14,rcut=30)
        for j in range(8,14):
            if abs(N[1][j][0]) > 4:
                xx.append(sorted([M[i],M[N[0][j]]]))
            if abs(N[1][j][1]) > 4:
                yy.append(sorted([M[i],M[N[0][j]]]))
            if abs(N[1][j][2]) > 4:
                zz.append(sorted([M[i],M[N[0][j]]]))
    # xx_hist[[same],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    def p_hist(xx):
        xx_hist = [0,0,0,0,0,0,0]
        for i in xx:
            if i[0]==i[1]:
                xx_hist[0] += 1
            if i == [0,1]:
                xx_hist[1] += 1
            if i == [0,2]:
                xx_hist[2] += 1
            if i == [0,3]:
                xx_hist[3] += 1
            if i == [1,2]:
                xx_hist[4] += 1
            if i == [1,3]:
                xx_hist[5] += 1
            if i == [2,3]:
                xx_hist[6] += 1
        for i in range(len(xx_hist)):
            xx_hist[i] = xx_hist[i]/2
        return xx_hist
    # op_hist[[same],[30],[60]]
    # green-blue, red-orange
    # 0 : blue, 1:red , 2:green, 3:orange
    def op_hist(xx):
        xx_hist = [0,0,0]
        for i in xx:
            if i[0]==i[1]:
                xx_hist[0] += 1
            if i == [1,3]:
                xx_hist[1] += 1
            if i == [0,2]:
                xx_hist[1] += 1
            if i == [0,1]:
                xx_hist[2] += 1
            if i == [0,3]:
                xx_hist[2] += 1
            if i == [1,2]:
                xx_hist[2] += 1
            if i == [2,3]:
                xx_hist[2] += 1
        for i in range(len(xx_hist)):
            xx_hist[i] = xx_hist[i]/2
        return xx_hist
    print 'pressure same twin different'
    #print '[same],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]'
    #print pxx,p_hist(xx)
    #print pyy,p_hist(yy)
    #print pzz,p_hist(zz)
    print pp
    print pxx, op_hist(xx)
    print pyy, op_hist(yy)
    print pzz, op_hist(zz)
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    for i in M:
        if i ==0:
            a1+=1
        if i ==1:
            a2+=1
        if i ==2:
            a3+=1
        if i ==3:
            a4+=1
    print a1,a2,a3,a4


#################################################
# xyz
#################################################
def run_single_xyz():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=readxyz.ReadCord(trajectory='C0_unit.xyz')
    L = np.array([[290., 290., 290.]])
    last = M.frames
    VW=M.cord_auto(['A'])
    x = [0]
    neighbor_no_gauss(VW, L, x,count=12, step=5e4,rcut=30)
############################################
def drift_remove(VW,L,step=5e4):
    from MD.analysis.drift_remove import eliminate_drift
    VW = eliminate_drift(VW,L)
    return VW
#################################################
# dcd
#################################################
def run_single():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    delta = 10
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    print_box_volume(L_cont,delta=delta)
    last = M.frames
    V_index = M.get_index(['V'])
    W_index = M.get_index(['W'])
    try:
        #V = util.pickle_load('V.pkl')
        Z = util.pickle_load('Z.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
    #    #V=M.cord_auto(['V'])
         Z=M.cord_auto(['Z'])
         VW=M.cord_auto(['V','W'])
         #VW,Z,W = drift_remove_all(VW,Z,VW,L)
    #    util.pickle_dump(V,'V.pkl')
         util.pickle_dump(Z,'Z.pkl')
         util.pickle_dump(VW,'VW.pkl')
    #Find the msd of the system
    #VW = drift_remove(VW,L)
    #write_xyz(VW)
    #x = range(0,last,delta)
    x = [0]
    L_last = L_cont[0][0]
    for i,j in enumerate(L_cont):
        if j[0] != L_last:
            L_last = j[0]
            if i - 5 > 0:
                x.append(i-5)
                #x.append(i)
    rotation_diffusion_split(Z,VW,L_cont,x)
    #x.append(last-5)
    #x.append(last)
    #a = 500
    #print x
    #for i in range(len(x)-1):
    #    a = min(a,x[i+1]-x[i])
    #for i in range(len(x)):
    #    x[i] = x[i]+a/2

    #del x[-1]
    #print x
    #print last
    #print a
    #del x[-1]
    #x = [last-10]
    #cubic_order_animate(VW,Z,L_cont,x)
    #linker_gauss_test(M,VW, Z,L, x, rcut=1.0, step=5e4)
    #linker_gauss(M,VW, Z,L, x, rcut=1.0, step=5e4)
    #cubic_order(VW,Z,L_cont,x)
    #cube_connected_rotate(VW,Z,M, L_cont)
    #delta = 10
    #x = range(0,VW.shape[0],delta)
    #x = [last-10,last-1]
    #x = range(last-10,last)
    #x = range(last-2,last)
    #x= [4]
    x = range(0,last,2)
    #cubic_order(VW,Z,L_cont,x)
    #frames = [25,last-5]
    #polymer length
    #center_end_connected(M,L_cont)
    #cube_connected_single_face(M, L_cont)
    #end_end_connected(M,L_cont)
    #end_end(M,L_cont)
    #cube_connected_angle(VW, Z, M, L_cont)
    #cube_angle(Z,VW,L_cont,delta=delta)
    #rotation_diffusion(Z,VW,L_cont)
    #######################
    ######################
    #msd(VW,L_cont[0],step=1)
    #rotations_area(Z,VW,L_cont)
    import MD.analysis.cubatic_order as c
    print VW.shape
    #assymetric_pressure(VW,L_cont,1710)
    #for i in range(0,200,10):
    #for i in range(101,200,5):
    #    try:
    #        assymetric_pressure(VW,L_cont,i)
    #    except:
    #        print 'error wrong frame number'
    #x = range(0,last-1,50)
    #polymer_interpenetrate(M,VW, Z,L_cont,x)
    #x = [0, 85,1700]
    #sp20
    #x = [200,350,400,435]
    #sp15
    #x = [40,140]
    #sp12
    #x = [230,325,1008,1445]
    #x = [0, 35,1000]
    #sp8
    #x = [80,586,793,1005]
    #sp4
    #x = [0,90,940]
    #x = [80,180]
    #x = [10,90]
    #x = 100,145
    #x = [0,30,230]
    #x = [15]
    #polymer_angle(M,VW,Z,L_cont,x)
    #polymer_azimuth(M,VW, Z,L_cont, x)
    #total_polymer_azimuth(M,L_cont,x)
    #polymer_rho(M,L_cont,x)
    #polymer_radial(M,L_cont,x)
    #total_polymer_radial(M,L_cont,x)
    #total_polymer_rho(M,L_cont,x)
    #draw_polymer(M,L_cont,x)
    #x = [985]
    #polymer_angle_hist(M,VW,Z,L_cont,x)
    #nn_distance(VW, L_cont, x, rcut=25.0)
    #x = range(0,last-10,10)
    print x
    #x = [184]
    #neighbor_gauss(VW, Z,L_cont, x, step=5e4,rcut=20)
    #neighbor_no_gauss(VW, L_cont, x,count=8, step=5e4,rcut=15.5)
    #neighbor_no_gauss(VW, L_cont, x,count=6 ,step=5e4,rcut=30)
    #neighbor_no_gauss(VW, L_cont, x,count=12 ,step=5e4,rcut=35)
    #neighbor_no_gauss(VW, L_cont, x,count=24 ,step=5e4,rcut=35)
    #x = range(0,last-5,5)
    #polymer_gauss(M,VW, Z,L_cont, x)
    #polymer_inter_gauss(M,VW, Z,L_cont, x,rcut=16)
    #ssDNA_gauss(M,VW, Z,L_cont, x)
    #c.gaussmap_cluster(VW,Z,L_cont,x,1,scale=12)
    #c.gaussmap(VW,Z,L_cont,x,1,scale=12)
    #t = 30
    #dt = 15
    #VW = Average(VW,L_cont,t,dt,save='star.xyz')
    #Z = Average(Z,L_cont,t,dt,save='Z.xyz')
    #c.gaussmap(np.array([VW]),np.array([Z]),L_cont,[0],1)
    #neighbor_no_gauss(np.array([VW]), L_cont, [0],count=12, step=5e4,rcut=25)
def cube_analysis():
    delta = 1
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    print_box_volume(L_cont,delta=delta)
    last = M.frames
    V_index = M.get_index(['V'])
    W_index = M.get_index(['W'])
    try:
        #V = util.pickle_load('V.pkl')
        Z = util.pickle_load('Z.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
    #    #V=M.cord_auto(['V'])
         Z=M.cord_auto(['Z'])
         VW=M.cord_auto(['V','W'])
         #VW,Z,W = drift_remove_all(VW,Z,VW,L)
    #    util.pickle_dump(V,'V.pkl')
         util.pickle_dump(Z,'Z.pkl')
         util.pickle_dump(VW,'VW.pkl')
    #Find the msd of the system
    #x = [0]
    #L_last = L_cont[0][0]
    #for i,j in enumerate(L_cont):
    #    if j[0] != L_last:
    #        L_last = j[0]
    #        if i - 5 > 0:
    #            x.append(i-5)
    #            #x.append(i)
    #x.append(last-5)
    x = [0]
    count = 0
    lcut = 5
    L_length = [[L_cont[0][0]]]
    L_last = L_cont[0]
    for i in range(1,len(L_cont),delta):
        if L_cont[i] != L_last :
            print i
            L_length[-1].append(count)
            L_length.append([L_cont[i][0]])
            count = 1
            L_last = L_cont[i]
        elif i >= len(L_cont)-delta:
            count+=1
            L_length[-1].append(count)
            break
        else:
            count+=1
    count = 0
    for i in L_length:
        if i[1] > lcut:
            x.append(i[1]+count)
            count+=i[1]
        else:
            print i[1]
            x[-1] += i[1]
            count+=i[1]

    print x
    #delta = 5
    #x = [0]
    #L_last = L_cont[0][0]
    #for i in range(0,len(L_cont),delta):
    #    if L_cont[i][0] != L_last or i > len(L_cont)-delta:
    #        L_last = L_cont[i][0]
    #        if i - 2 > 0:
    #            x.append(i-1)
    #print x
    #rotation_diffusion_split(Z,VW,L_cont,x)
    #delta = 3
    #x = range(0,last,delta)
    ##cubic_order(VW,Z,L_cont,x)
    #x = [last-10]
    #xx = []
    #for i in x[1:]:
    #    xx.append(i-10)
    #print x
    #print xx
    #xx = [40]
    #xx = [150,237,290,344,390]
    xx=[150]
    #t = 150
    #dt = 5
    #VW = Radial_Average(VW,VW[t][1010],L_cont,t,dt,save='radial_star.xyz')
    neighbor_no_gauss(VW, L_cont, xx,count=8, step=5e4,rcut=25.5)
    import MD.analysis.cubatic_order as c
    #t = 150
    #dt = 5
    #VW = Average_simple(VW,L_cont,t,dt,save='star.xyz')
    #delta =5
    #x = range(0,last,delta)
    #c.gaussmap(VW,Z,L_cont,x,1,scale=12)
##For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            try:
#                run_single()
#            except:
#                pass
#            print '##\nfinished with:',f
#            os.chdir('../')
#For single directories
if __name__ == '__main__':
    #run_debug()
    #run_all()
    #run_single()
    #run_single()
    #krun_simple()
    #run_binary()
    cube_analysis()

