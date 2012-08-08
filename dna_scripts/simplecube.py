
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
from MD.analysis.particle_distance import particle_distance
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import second_nearest_neighbors_index
from MD.analysis.nearest_neighbor import min_distance
from MD.plot.histogram import histogram
from MD.plot.histogram import histogram_normal
from MD.plot.histogram import histogram_reg
from MD.analysis.rotation import single_rotation
############ 
# Run from inside a folder
#########
reload(MD)
reload(pyplot)

########################################3
def dump_xyz(V):
    out = open('simple.xyz','w')
    print V.shape
    for i in range(V.shape[0]):
        out.write(('%i\n\n')%(V.shape[1]))
        for j in range(V.shape[1]):
            out.write(('V %.2f %.2f %.2f\n')%(V[i][j][0],V[i][j][1],V[i][j][2]))
############################################
## msd
############################################
def msd(VW,L,step=1):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x,msd=msd(VW,L,step=step)
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save='MSDtime')
    util.pickle_dump(msd,'msd.pkl')
############################################
## returns the arrays with the drift removed
############################################
def drift_remove(VW,L,index):
    from MD.analysis.drift_remove import eliminate_drift
    VW = eliminate_drift(VW,L,index)
    return VW
    #Find the msd of the system
############################################
## returns the arrays with the drift removed
############################################
def drift_remove_all(VW,V,W,L):
    from MD.analysis.drift_remove import eliminate_drift
    W = eliminate_drift(W,L)
    V = eliminate_drift(V,L)
    VW = eliminate_drift(VW,L)
    return VW,V,W
    #Find the msd of the system
############################################
## 	g(r)
############################################
def distance_distribution(V,L,n_frames=15,n_start=1,end=10):
    #finds end to end distances and plots a histogram with that data
    def hist(M1,M2,L,save_name,label,end=1):
        distance=particle_distance(M1,M2,L)
        hist_s,xs,max_hist=histogram_normal(distance,bins=60,end=end)
        #normalize the function with respect to an ideal gas
        delta_r = xs[1] - xs[0]
        for i in range(len(hist_s)):
            r = (xs[i])
            hist_s[i] /= (4*math.pi*delta_r*r*(M1.shape[1]/L[0]**3))
        return xs,hist_s
    ########################################
    #Plot Data for AA,AB,BB distances
    save_name='_time_'+'%i'%(n_start)
    start=n_start
    finish=n_start+n_frames
    max_hist=0
    BB=hist(V[start:finish],V[start:finish],L,save_name,'B-B',end=end)
    #q is the first peak in the plot or the lattice spacing
    q=13
    pyplot.plot_sf_sc(BB[0][1:],BB[1][1:],q,xlabel='s',ylabel='g(s)',save='nnplot_'+save_name)
    out = open('dist.txt','w')
    out.write('\n')
    for i in BB[1]:
        out.write(('%.2f\n')%(i))
    out.close()

############################################
# lifetime of connections
########################################
def find_lifetime(M,L,steps,temp,step_range=30,delta=4,rcut=1.0,step=5e4):
    import MD.analysis.lifetime as life
    try:
        C = util.pickle_load('C.pkl')
        G = util.pickle_load('G.pkl')
    except:
        C=M.cord_auto(['C'])
        G=M.cord_auto(['G'])
        util.pickle_dump(C,'C.pkl')
        util.pickle_dump(G,'G.pkl')
    #The total number of frames we are going to look at
    for i in steps:
        print 'Frame',i,'Temp',temp[i]
        x=np.arange(i,step_range+i,delta)
        #Find he number of connections at specific points x
        remain = life.lifetime(C[x],G[x],L)
        print remain
        pyplot.plot(x,remain,xlabel='Time',
                ylabel='remaining connections', save='lifetime'+temp[i])
        plt.close()
######################
# FInd the solid particles in a simulation
#####################
def solid_particles_vmd(VW,V,W,L,var,skip=15,step=5e4,rcut=False):
    import MD.dna_scripts.solidparticle as sp
    reload(sp)
    if os.path.exists('sccrystals.xyz') == False:
        sp.solid_particles_sc(VW,V,W,L,skip=skip,step=5e4,rcut=False)
    sp.crystal_size_cube(L,VW.shape[0]/skip,VW.shape[1])
######################
# FInd the Angle between cubes
#####################
def cube_angle(M,L):
    Z = M.cord_auto(['Z'])
    VW = M.cord_auto(['V','W'])
    ### We need to know which z vector to use
    ### We can find the minimum distance to Z bead from M bead and use that one
    ### to identify the side of the cube that the ssDNA is attatched. 
    def min_distance(Z1,Z2,L):
        num=0
        side1=0
        side2=0
        smallest=50
        for i in range(Z1.shape[0]):
            for j in range(Z2.shape[0]):
                num = points.dist(Z1[i],Z2[j],L)[0]
                if smallest > num:
                    smallest = num
                    side1 = i
                    side2 = j
        return side1,side2

    rcut=min_particle_distance(VW[-1:],VW[-1:],L)+6
    #loop through and azimithal distance
    theta=[]
    for k in range(VW.shape[0]-1,VW.shape[0]):
        print k
        for i in range(VW.shape[1]):
            NN=nearest_neighbors_index(VW[k],i,L,rcut=rcut)
            for j in NN[0]:
                #project vector from center to end onto z vector
                side1,side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                v1 = points.vector(VW[k][i],Z[k][side1+i*6],L)
                v2 = points.vector(VW[k][j],Z[k][side2+j*6],L)
                theta.append(points.angle_between(v1,v2))
    hist_theta,xtheta,max_theta = histogram(theta,50)
    pyplot.plot_bar(xtheta,hist_theta,
            xlabel='theta radaians',ylabel='amount',save='nntheta')
    return xtheta, hist_theta
######################
# FInd the Angle between cubes
#####################
def cube_second_angle(M,L):
    Z = M.cord_auto(['Z'])
    VW = M.cord_auto(['V','W'])
    ### We need to know which z vector to use
    ### We can find the minimum distance to Z bead from M bead and use that one
    ### to identify the side of the cube that the ssDNA is attatched. 
    def min_distance(Z1,Z2,L):
        num=0
        side1=0
        side2=0
        smallest=50
        for i in range(Z1.shape[0]):
            for j in range(Z2.shape[0]):
                num = points.dist(Z1[i],Z2[j],L)[0]
                if smallest > num:
                    smallest = num
                    side1 = i
                    side2 = j
        return side1,side2

    rcut=min_particle_distance(VW[-1:],VW[-1:],L)+6
    #loop through and azimithal distance
    theta=[]
    for k in range(VW.shape[0]-1,VW.shape[0]):
        print k
        for i in range(VW.shape[1]):
            NN  = second_nearest_neighbors_index(VW[k],i,L,inrcut=rcut*3/2,outrcut=2*rcut+4)
            for j in NN[0]:
                #project vector from center to end onto z vector
                side1,side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                v1 = points.vector(VW[k][i],Z[k][side1+i*6],L)
                v2 = points.vector(VW[k][j],Z[k][side2+j*6],L)
                theta.append(points.angle_between(v1,v2))
    hist_theta,xtheta,max_theta = histogram(theta,50)
    pyplot.plot_bar(xtheta,hist_theta,
            xlabel='theta radaians',ylabel='amount',save='second_nntheta')
    return xtheta, hist_theta
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

## \brief find the anlge between cubes we know have connections
#
# \returns x,y values of histogram (theat, # at theta)
#
# \param 
# \param VW
# \param L
# \param ndna number of dna per cube
# \param frame frame to look at
def cube_rotated_connected(M, L, ndna,frame=-1):
    Z = M.cord_auto(['Z'])
    VW = M.cord_auto(['V','W'])
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    try:
        connections=util.pickle_load('concube.pkl')
    except:
        C=M.cord_auto(['C'])
        G=M.cord_auto(['G'])
        connections = con.connections(np.array([C[frame]]),np.array([G[frame]]),L)
        util.pickle_dump(connections,'concube.pkl')
    #find a graph of the particles
    import MD.analysis.graph as graph
    import networkx as nx
    networks, num_networks, deg, neighbors, num_nieghbors, gr = graph.grapher(connections,VW.shape[1],ndna)
    ##################################
    #loop through and azimithal distance
    def make_plot(VW,Z,gr,equal=1,k=1):
        def find_side(side):
            if side%6 + 2 > 6:
                s2 = side - 2
                s3 = side - 2
            else:
                if side%6 - 2 < 0:
                    s2 = side + 2
                    s3 = side + 2
                else:
                    s2 = side + 2
                    s3 = side - 2
            return s2,s3
        dh = 25
        theta=[]
        count=0
        for i in range(VW.shape[1]):
            #get a list of the particles which are connected
            NN = gr.adjacency_list()
            for j in NN[i]:
                #project vector from center to end onto z vector
                if gr.number_of_edges(i,j) == equal:
                    side1, side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                    C = VW[k][i]
                    s2, s3 = find_side(side1)
                    v1 = points.vector( points.apply_vboundary(VW[k][i]-C,L),
                                        points.apply_vboundary(Z[k][side1+i*6]-C,L), L)
                    z1 = points.vector( points.apply_vboundary(VW[k][i]-C,L),
                                        points.apply_vboundary(Z[k][s2+i*6]-C,L), L)
                    v2 = points.vector( points.apply_vboundary(VW[k][j]-C,L),
                                        points.apply_vboundary(Z[k][side2+i*6]-C,L), L)
                    #get points at 0
                    proj = points.vector_plane_projection(v2,v1)
                    print v1
                    print z1
                    print v2
                    print proj
                    theta.append(points.angle_between(z1, proj))
                    count += 1
        try:
            hist_theta, xtheta = histogram_reg(theta,dh)
            pyplot.plot_bar(xtheta,hist_theta,
                    xlabel='theta radaians',ylabel='amount',
                    save='rotate%i_count%i'%(equal,count), width=0.05)
        except:
            pass
        plt.close()
    #####################################
    k = frame
    #plot which takes into account the number of connections
    make_plot(VW,Z,gr,equal=1,k=frame)
    make_plot(VW,Z,gr,equal=2,k=frame)
    make_plot(VW,Z,gr,equal=3,k=frame)

    def find_side(side):
        if side%6 + 2 > 6:
            s2 = side - 2
            s3 = side - 2
        else:
            if side%6 - 2 < 0:
                s2 = side + 2
                s3 = side + 2
            else:
                s2 = side + 2
                s3 = side - 2
        return s2,s3
    #plot which takes into account the number of connections
    theta=[]
    count=0
    dh = 25
    for i in range(VW.shape[1]):
        #get a list of the particles which are connected
        NN = gr.adjacency_list()
        for j in NN[i]:
            #project vector from center to end onto z vector
            if gr.number_of_edges(i,j) > 3:
                side1, side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                C = VW[k][i]
                s2, s3 = find_side(side1)
                v1 = points.vector( points.apply_vboundary(VW[k][i]-C,L),
                                    points.apply_vboundary(Z[k][side1+i*6]-C,L), L)
                z1 = points.vector( points.apply_vboundary(VW[k][i]-C,L),
                                    points.apply_vboundary(Z[k][s2+i*6]-C,L), L)
                v2 = points.vector( points.apply_vboundary(VW[k][j]-C,L),
                                    points.apply_vboundary(Z[k][side2+i*6]-C,L), L)
                #get points at 0
                proj = points.vector_plane_projection(v2,v1)
                theta.append(points.angle_between(z1, proj))
                count += 1
    hist_theta,xtheta = histogram_reg(theta,dh)
    print theta
    print hist_theta
    pyplot.plot_bar(xtheta,hist_theta,
            xlabel='theta radaians', ylabel='amount',
            save='rotate_g3_count%i'%count,width=0.05)
    plt.close()
    return xtheta, hist_theta
## \brief find the anlge between cubes we know have connections
#       but look at the second connection
#
# \returns x,y values of histogram (theat, # at theta)
#
# \param 
# \param VW
# \param L
# \param ndna number of dna per cube
# \param frame frame to look at
def second_cube_angle_connected(M, L, ndna,frame=-1):
    Z = M.cord_auto(['Z'])
    VW = M.cord_auto(['V','W'])
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    try:
        connections=util.pickle_load('concube.pkl')
    except:
        C=M.cord_auto(['C'])
        G=M.cord_auto(['G'])
        connections = con.connections(np.array([C[frame]]),np.array([G[frame]]),L)
        util.pickle_dump(connections,'concube.pkl')
    #find a graph of the particles
    ##################################
    #loop through and azimithal distance
    def make_plot(connections, VW, Z, ndna, k=1,f1=1,f2=1, greater = False):
        import MD.analysis.graph as graph
        reload(graph)
        dh = 25
        theta=[]
        count=0
        #get a list of the particles which are connected
        NN = graph.second_grapher(connections, VW.shape[1], ndna,
                f1cut = f1, f2cut = f2, greater = greater)
        if NN != []:
            for i in range(VW.shape[1]):
                for j in NN[i]:
                    #project vector from center to end onto z vector
                    side1,side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                    v1 = points.vector(VW[k][i],Z[k][side1+i*6],L)
                    v2 = points.vector(VW[k][j],Z[k][side2+j*6],L)
                    theta.append(points.angle_between(v1,v2))
                    count += 1
            hist_theta, xtheta = histogram_reg(theta,dh)
            pyplot.plot_bar(xtheta,hist_theta,
                    xlabel='theta radaians',ylabel='amount',
                    save='theta%i_%i_count%i'%(f1, f2, count), width=0.05)
            plt.close()
        else:
            print NN
    #####################################
    k = frame
    #plot which takes into account the number of connections
    make_plot(connections[0], VW, Z, ndna, k = frame, f1 = 1, f2 = 1)
    make_plot(connections[0], VW, Z, ndna, k = frame, f1 = 2, f2 = 2)
    make_plot(connections[0], VW, Z, ndna, k = frame, f1 = 3, f2 = 3)
    make_plot(connections[0], VW, Z, ndna, k = frame, f1 = 4, f2 = 3,
            greater = True)

## \neighbors not connected
#
# \returns x,y values of histogram (theat, # at theta)
#
# \param 
# \param M
# \param L
# \param ndna number of dna per cube
# \param frame frame to look at
def not_connected_nn(M,VW, L, ndna, n_finish=1, n_start=0, delta=1, rcut = False):
    C=M.cord_auto(['C'])
    G=M.cord_auto(['G'])
    #The total number of frames we are going to look at
    n_frames = (n_finish-n_start)/delta
    x=np.arange(n_start,n_finish,delta)
    import MD.analysis.connections as con
    if rcut == False:
        rcut=min_particle_distance(VW[-1:],VW[-1:],L)+3
    try:
        connections=util.pickle_load('con.pkl')
    except:
        C=M.cord_auto(['C'])
        G=M.cord_auto(['G'])
        connections = con.connections(np.array([C[x]]),np.array([G[x]]),L)
        util.pickle_dump(connections,'con.pkl')
    #find a graph of the particles
    import MD.analysis.graph as graph
    import networkx as nx
    ##################################
    #loop through and azimithal distance
    pin = 0
    counter = []
    net, num, deg, NC, num_n, gr = graph.grapher(connections,VW.shape[1],ndna)
    for k in x:
        print k
        count = 0
        for i in range(len(NC[pin])):
            NN=nearest_neighbors_index(VW[k],i,L,rcut=rcut)[0]
            print NN
            print NC[pin][i]
            count += len(NN)-len(NC[pin][i])
        counter.append(count)
        pin+=1
    pyplot.plot(x,np.array(counter),xlabel='Time',
            ylabel='nn not connected', save='nnnotcon')
    return x, counter
## \brief find the number of cubes which are connected from dna at one or more 
#       faces
#
# \returns x,y1,y2,y3,y4,y5 values 
#
# \param 
# \param M
# \param L
# \param ndna number of dna per cube
# \param n_nodes number of Nanoparticles in the system
def cube_connected_face(M, L, ndna, n_nodes, n_start = 0, n_finish = 1, delta = 30):
    x=np.arange(n_start,n_finish,delta)
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    reload(graph)
    try:
        connections = util.pickle_load('con.pkl')
        delta = n_finish/len(connections)
        x=np.arange(n_start,n_finish,delta)
        x=np.arange(len(connections))
    except:
        G=M.cord_auto(['G'])
        C=M.cord_auto(['C'])
        connections = con.connections(C[x],G[x],L,rcut=1.0)
        util.pickle_dump(connections,'con.pkl')
        util.pickle_dump(C,'C.pkl')
        util.pickle_dump(G,'G.pkl')
    counter = graph.sides(connections,n_nodes,ndna,cut=5)
    counter1 = graph.sides(connections,n_nodes,ndna,cut=4)
    counter2 = graph.sides(connections,n_nodes,ndna,cut=3)
    counter3 = graph.sides(connections,n_nodes,ndna,cut=2)
    counter4 = graph.sides(connections,n_nodes,ndna,cut=6)
    pyplot.plot(x,counter, xlabel='time',ylabel='amount',
            save='faces_1')
    pyplot.plot(x,counter1, xlabel='time',ylabel='amount',
            save='faces_2')
    pyplot.plot(x,counter2, xlabel='time',ylabel='amount',
            save='faces_3')
    pyplot.plot(x,counter3, xlabel='time',ylabel='amount',
            save='faces_4')
    pyplot.plot(x,counter4, xlabel='time',ylabel='amount',
            save='faces_0')
    return x, counter, counter1, counter2, counter3
## \brief find the number of cubes connected to a single face on another cube
#
# \returns x,y1,y2,y3,y4,y5 values 
#
# \param 
# \param M
# \param L
# \param ndna number of dna per cube
# \param n_nodes number of Nanoparticles in the system
def cube_connected_single_face(M, L, ndna, n_nodes, n_start = 0, n_finish = 1, delta = 30):
    x=np.arange(n_start,n_finish,delta)
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    reload(graph)
    try:
        connections = util.pickle_load('con.pkl')
        delta = n_finish/len(connections)
        x=np.arange(len(connections))
        x=np.arange(len(connections))
    except:
        G=M.cord_auto(['G'])
        C=M.cord_auto(['C'])
        connections = con.connections(C[x],G[x],L,rcut=1.0)
        util.pickle_dump(connections,'con.pkl')
    counter0 = graph.sides_mod(connections,n_nodes,ndna,cut=0)
    counter1 = graph.sides_mod(connections,n_nodes,ndna,cut=1)
    counter2 = graph.sides_mod(connections,n_nodes,ndna,cut=2)
    counter3 = graph.sides_mod(connections,n_nodes,ndna,cut=3)
    counter4 = graph.sides_mod(connections,n_nodes,ndna,cut=4)
    counter5 = graph.sides_mod(connections,n_nodes,ndna,cut=5)
    counter5 = graph.sides_mod(connections,n_nodes,ndna,cut=6)
    pyplot.plot(x,counter0, xlabel='time',ylabel='amount',
            save='confaces_0')
    pyplot.plot(x,counter1, xlabel='time',ylabel='amount',
            save='confaces_1')
    pyplot.plot(x,counter2, xlabel='time',ylabel='amount',
            save='confaces_2')
    pyplot.plot(x,counter3, xlabel='time',ylabel='amount',
            save='confaces_3')
    pyplot.plot(x,counter4, xlabel='time',ylabel='amount',
            save='confaces_4')
    pyplot.plot(x,counter5, xlabel='time',ylabel='amount',
            save='confaces_5')
    return x, counter0, counter1, counter2, counter3
######################
# FInd the Angle between cube normal and c-axis
# caxis is the vector of the caxis
#####################
def caxis_angle(M, VW, L, caxis, frame=[-1], sides=6):
    Z = M.cord_auto(['Z'])
    #loop through and azimithal distance
    theta_1=[]
    theta_2=[]
    for k in frame:
        print 'frame',k
        for i in range(VW.shape[1]):
            min_theta=100
            side = 0
            for j in range(sides):
                #project vector from center to end onto z vector
                v1 = points.vector(VW[k][i],Z[k][j+i*sides],L)
                theta = points.angle_between(v1,caxis)
                if min_theta>theta:
                    min_theta = theta
                    side = j
            if side>4:
                v2 = points.vector(VW[k][i],Z[k][side-2+i*sides],L)
                theta2 = points.angle_between(v2,caxis)
            if side<2:
                v2 = points.vector(VW[k][i],Z[k][side+2+i*sides],L)
                theta2 = points.angle_between(v2,caxis)
            if side <= 4 and side >= 2:
                v2 = points.vector(VW[k][i],Z[k][1+2+i*sides],L)
                theta2 = points.angle_between(v2,caxis)
            theta_1.append(min_theta)
            theta_2.append(theta2)
    #normalize
    hist_theta, xtheta, max_theta = histogram_normal(theta_1,8)
    hist_theta2, xtheta2, max_theta2 = histogram_normal(theta_2,8)
    import math
    print xtheta
    for i in range(hist_theta.shape[0]):
        print hist_theta[i]
        hist_theta[i] /= math.sin((xtheta[i]+xtheta[i+1])/2.0)
        hist_theta2[i] /= math.sin((xtheta2[i]+xtheta2[i+1])/2.0)
    pyplot.plot_bar(xtheta[1:],hist_theta,
            xlabel='theta radaians',ylabel='amount',save='caxis_theta')
    plt.close()
    pyplot.plot_bar(xtheta2[1:],hist_theta2,
            xlabel='theta radaians',ylabel='amount',save='caxis_theta2')
    plt.close()
    return xtheta, hist_theta,  xtheta2, hist_theta2
######################
# FInd the rotational motion of cube
#####################
def rotations(M,VW,L,chord):
    x=range(0,VW.shape[1],10)
    Z = M.cord_auto([chord])
    for i in x:
        theta, phi = single_rotation(VW,Z,L,index=i)
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
#####################################################
#find number of a specified pair of particles that are within r_cut of each
#####################################################
def find_networks(M, VW, L, var, n_finish=1, n_start=0, delta=30, rcut=1.0):
    C=M.cord_auto(['C'])
    ndna = (C.shape[1]*2) / VW.shape[1]
    #The total number of frames we are going to look at
    x=np.arange(n_start,n_finish,delta)
    print len(x)
    #Find the number of connections at specific points x
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    try:
        connections = util.pickle_load('con.pkl')
        delta = n_finish/len(connections)
        x=np.arange(n_start,n_finish,delta)
    except:
        G=M.cord_auto(['G'])
        connections = con.connections(C[x],G[x],L,rcut=rcut)
        util.pickle_dump(connections,'con.pkl')
    #plot total number of connections
    print connections
    num_connections=con.num_connections(connections,VW.shape[1])
    pyplot.plot(x,num_connections,xlabel='Time',
            ylabel='hyrbid. density', save='connections')
    plt.close()
    #get the info
    networks, num_networks, deg, neighbors, num_n, gr = graph.grapher(connections,VW.shape[1],ndna)
    #plot the number of neighbors at each timesteps
    pyplot.plot(x,num_n,xlabel='Time',
            ylabel='num neighbors', save='neighbors')
    plt.close()
    print 'making plot'
    label = ['networks','1','2','3','4','5','6','7','8','9']
    pyplot.plot_multi(x,deg,label,xlabel='time',ylabel='number',save='con_net')
    return x, networks, connections
#####################################################
# Find Static Structure Facto
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor(VW, L, n_frames=10,
        n_start=1,l=10,filters=0.05,save='sf'):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    stmp,qx = sf.sfactor(VW[n_start:n_start+n_frames],L=L,l=l)
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
    #Plot the graph
    xlabel = '$|\\vec{q}$ $|$'
    ylabel = '$S(\\vec{q}$ $)$ '
    pyplot.plot(Q, S, linestyle='', marker='x',
            xlabel=xlabel, ylabel=ylabel, save=save)
    util.pickle_dump(Q,'sfx'+str(n_start)+'.pkl')
    util.pickle_dump(S,'sfy'+str(n_start)+'.pkl')
    #if len(QQ)>500:
    #    bars,height = sf.sf_max(QQ,qpeaks=8)
    #    db_x,db_y = sf.dbfactor(msdr=5)
    #    pyplot.plot_bar(bars, height, xlabel=xlabel, ylabel=ylabel,
    #            save=save)
    #    pyplot.plot(db_x, db_y, linestyle='--', xlabel=xlabel, ylabel=ylabel,
    #            save=save)
    #    #return q, s, bars, height, db_x, db_y
    #    return Q, S, primitive_vectors
    #filter out the noise
    Q, S, primitive_vectors = sf.sf_filter(Q, S, primitive_vectors,
            filters=filters, save = save)
    return Q, S, primitive_vectors
#####################################################
# Look at mixing of species of spheres
#####################################################
def species_mixing(V,W,L,n_finish=1,n_start=0,delta=100,rcut=False):
    #The total number of frames we are going to look at
    if rcut == False:
        rcut=min_particle_distance(V[-1:],W[-1:],L)+6
    print rcut
    #rcut = 19
    rcut = 20
    n_frames = (n_finish-n_start)/delta
    x = np.arange(n_start,n_finish,delta)
    s = []
    o = []
    for k in x:
        Same = 0
        Other = 0
        print k
        for i in range(V[k].shape[0]):
            Same += len(nearest_neighbors_point(V[k],V[k][i],L,rcut)[0])
            Same += len(nearest_neighbors_point(W[k],W[k][i],L,rcut)[0])
            Other += len(nearest_neighbors_point(W[k],V[k][i],L,rcut)[0])
            Other += len(nearest_neighbors_point(V[k],W[k][i],L,rcut)[0])
        s.append(float(Same)/(W[k].shape[0]+V[k].shape[0]))
        o.append(float(Other)/(V[k].shape[0]+W[k].shape[0]))
    print s
    print o
    util.pickle_dump([x,s,o],'mix.pkl')
    pyplot.plot2(x, s, x, o, label1='Same',
            label2='other', save='mixing_avg', showleg=True)
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
## \brief read mylog file and makes a plot of the energy
#
# \returns x coordinates of time step
# \returns y whatever value of the row you are looking for
#
# \param row the row you would like to grab data from 
def mylog_volume(M, num= 128, row = 2):
    aasdf
# end to end distance
###################3
def end_end(V,W,M,L,var):
    T=M.cord_auto(['T'])
    K=M.cord_auto(['K'])
    V=M.cord_auto(['V'])
    W=M.cord_auto(['W'])
    ndna = T.shape[1]/V.shape[1]
    #K-V
    #W-T
    VK=[]
    WT=[]
    last = V.shape[0]
    #find center-ssdna particle distance
    for k in range(last-5,last):
        for i in range(V.shape[1]):
            for j in range(ndna):
                VK.append(points.dist(V[k][i],K[k][j+i*ndna],L)[0])
                WT.append(points.dist(W[k][i],T[k][j+i*ndna],L)[0])
    hist_K,xK,max_k=histogram(VK,50)
    hist_T,xT,max_T=histogram(WT,50)

    pyplot.plot_bar(xK,hist_K,save='VK_distance')
    plt.close()
    pyplot.plot_bar(xT,hist_T,save='WT_distance')
    plt.close()

    #find the length of the polymer

    KT=M.cord_auto(['K','T'])
    S=M.cord_auto(['M'])
    edge=[]
    for k in range(last-5,last):
        for i in range(S.shape[1]):
            edge.append(points.dist(KT[k][i],S[k][i],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,50)

    pyplot.plot_bar(xedge,hist_edge,save='polymer_length')
    plt.close()
## find the z directional length of polymer
## away from the surface of the cube
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
# What to run
#################################################
def run_test():
    #setup
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = 38.71
    Lx = L
    Ly = L
    Lz = L
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    last = M.frames
    try:
        V = util.pickle_load('V.pkl')
        Z = util.pickle_load('Z.pkl')
    except:
        V=M.cord_auto(['V'])
        V=M.cord_auto(['Z'])
        dump_xyz(V)
        dump_xyz(Z)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(Z,'V.pkl')
    delta = last/5
    x = range(last-delta,last,delta)
    #x = [0]
    #C = L_cont[0][0]
    #use when changing the box size
    #for i in range(1,last-2):
    #    if L_cont[i][0] != C:
    #        x.append(i)
    #        C = L_cont[i][0]
    #        print i
    #        print C
    #cubics(M,VW,L,var['ndna'])
    #plt.close()
    #step,temp = mylog(row = 2)
    #plt.close()
    #rotations(M,V,L,'Z')
    #plt.close()
    #msd(V,L)
    #plt.close()
    #mylog_volume(M, num=V.shape[1])
    #plt.close()
    #solid_particles_vmd(VW,V,W,L,var,skip=100)
    #plt.close()
    #s = [x[-1]]
    for i in x:
        #print "finding s(q)"
        #if i >150:
        #    structure_factor(V, L_cont[i], n_start=i+10,save='sf'+str(i))
        #    plt.close()
        #    print "finding g(s)"
        #    distance_distribution(V,L_cont[i],n_start=i+10)
        #    plt.close()
        #else:
        #structure_factor(V, L_cont[i], n_start=i+10,save='sf'+str(i),n_frames =5,l=20)
        structure_factor(V, L, n_start=i,save='sf'+str(i),n_frames = 10,l=10)
        #plt.close()
        print "finding g(s)"
        distance_distribution(V,L_cont[i],n_start=i)
        plt.close()
    #end_end(V,W,M,L,var)
    #cube_second_angle(M,L)
    #not_connected_nn(M,VW, L, var['ndna'], n_finish=last, n_start=0, delta=60)
    #second_cube_angle_connected(M, L, var['ndna'])
    #cube_rotated_connected(M, L, var['ndna'],frame = -1)
#For multiple directories
if __name__ == '__main__':
    for f in sorted(os.listdir("./")):
        if os.path.isdir(f):
            os.chdir(f)
            try:
                run_test()
                print '##\nfinished with:',f
            except:
                print 'failed in this directory'
            os.chdir('../')
#For single directories
if __name__ == '__main__':
    #run_debug()
    #run_all()
    #run_single()
    run_test()

