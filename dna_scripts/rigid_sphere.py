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
#########################################
# What you want to do in those directories
########################################
def end_end(V,W,M,L,var):
    T=M.cord_auto(['T'])
    K=M.cord_auto(['K'])
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
############################################
## msd
############################################
def msd(VW,L,step=1):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x,msd=msd(VW,L,step=step)
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save='MSDtime')
############################################
## returns the arrays with the drift removed
############################################
def drift_remove(VW,L):
    from MD.analysis.drift_remove import eliminate_drift
    VW = eliminate_drift(VW,L)
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
def distance_distribution(V,W,L,n_frames=15,n_start=1,end=10):
    #finds end to end distances and plots a histogram with that data
    def hist(M1,M2,save_name,label,end=1):
        distance=particle_distance(M1,M2,L)
        hist_s,xs,max_hist=histogram_normal(distance,bins=50,end=end)
        return xs,hist_s
    ########################################
    #Plot Data for AA,AB,BB distances
    save_name='_time_'+'%i'%(n_start)
    start=n_start
    finish=n_start+n_frames
    max_hist=0
    AB=hist(V[start:finish],W[start:finish],save_name,'A-B',end=end)
    AA=hist(W[start:finish],W[start:finish],save_name,'A-A',end=end)
    BB=hist(V[start:finish],V[start:finish],save_name,'B-B',end=end)
    pyplot.plot3(AA[0],AA[1],AB[0],AB[1],BB[0],BB[1],'s','g(s)',
            'A-A','A-B','B-B',save='nnplot_'+save_name,showleg=True)
############################################
# lifetime of connections
########################################
def find_lifetime(M,L,n_finish=30,n_start=1,delta=4,rcut=1.0,step=5e4):
    C=M.cord_auto(['C'])
    G=M.cord_auto(['G'])
    #The total number of frames we are going to look at
    n_frames = (n_finish-n_start)/delta
    x=np.arange(n_start,n_finish,delta)
    print x
    #Find he number of connections at specific points x
    import MD.analysis.lifetime as life
    remain = life.lifetime(C[x],G[x],L)
    print remain
    pyplot.plot(x,remain,xlabel='Time',
            ylabel='remaining connections', save='lifetime')
    return x,remain
######################
# FInd the solid particles in a simulation
#####################
def solid_particles_vmd(VW,V,W,L,var,skip=15,step=5e4,rcut=False):
    import MD.dna_scripts.solidparticle as sp
    reload(sp)
    if os.path.exists('bcccrystals.xyz') == False:
        sp.solid_particles(VW,V,W,L,skip=skip,step=5e4,rcut=False)
    #sp.crystal_size(L,VW.shape[0]/skip,VW.shape[1])
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
## \brief find the anlge between cubes we know have connections
#
# \returns x,y values of histogram (theat, # at theta)
#
# \param 
# \param VW
# \param L
# \param ndna number of dna per cube
# \param frame frame to look at
def cube_angle_connected(M, L, ndna,frame=-1):
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
        dh = 25
        theta=[]
        count=0
        for i in range(VW.shape[1]):
            #get a list of the particles which are connected
            NN = gr.adjacency_list()
            for j in NN[i]:
                #project vector from center to end onto z vector
                if gr.number_of_edges(i,j) == equal:
                    side1,side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                    v1 = points.vector(VW[k][i],Z[k][side1+i*6],L)
                    v2 = points.vector(VW[k][j],Z[k][side2+j*6],L)
                    theta.append(points.angle_between(v1,v2))
                    count += 1
        try:
            hist_theta, xtheta = histogram_reg(theta,dh)
            pyplot.plot_bar(xtheta,hist_theta,
                    xlabel='theta radaians',ylabel='amount',
                    save='theta%i_count%i'%(equal,count), width=0.05)
        except:
            pass
        plt.close()
    #####################################
    k = frame
    #plot which takes into account the number of connections
    make_plot(VW,Z,gr,equal=1,k=frame)
    make_plot(VW,Z,gr,equal=2,k=frame)
    make_plot(VW,Z,gr,equal=3,k=frame)

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
                side1,side2 = min_distance(Z[k][i*6:(i+1)*6],Z[k][j*6:(j+1)*6],L)
                v1 = points.vector(VW[k][i],Z[k][side1+i*6],L)
                v2 = points.vector(VW[k][j],Z[k][side2+j*6],L)
                theta.append(points.angle_between(v1,v2))
                count += 1
    hist_theta,xtheta = histogram_reg(theta,dh)
    print theta
    print hist_theta
    pyplot.plot_bar(xtheta,hist_theta,
            xlabel='theta radaians', ylabel='amount',
            save='theta_g3_count%i'%count,width=0.05)
    plt.close()
    return xtheta, hist_theta
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
        x = range(0,n_finish,n_finish/len(connections))
    except:
        G=M.cord_auto(['G'])
        C=M.cord_auto(['C'])
        connections = con.connections(C[x],G[x],L,rcut=1.0)
        util.pickle_dump(connections,'con.pkl')
    counter = graph.sides(connections,n_nodes,ndna,cut=5)
    counter1 = graph.sides(connections,n_nodes,ndna,cut=4)
    counter2 = graph.sides(connections,n_nodes,ndna,cut=3)
    counter3 = graph.sides(connections,n_nodes,ndna,cut=2)
    pyplot.plot(x,counter, xlabel='time',ylabel='amount',
            save='faces_1')
    pyplot.plot(x,counter1, xlabel='time',ylabel='amount',
            save='faces_2')
    pyplot.plot(x,counter2, xlabel='time',ylabel='amount',
            save='faces_3')
    pyplot.plot(x,counter3, xlabel='time',ylabel='amount',
            save='faces_4')
    return x, counter, counter1, counter2, counter3
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
######################
# FInd the rotational motion of sphere
#####################
def sphere_rotations(M,VW,L,chord):
    x=range(0,VW.shape[1],10)
    Z = M.cord_auto([chord])
    ZZ = np.zeros((Z.shape[0],60*VW.shape[1],3))

    print Z.shape[0]*Z.shape[1]
    print ZZ.shape[0]*ZZ.shape[1]
    print Z.shape[1]
    print ZZ.shape[1]
    print 96+60*VW.shape[1]
    print 96*VW.shape[1]
    print 60*VW.shape[1]
    for i in range(VW.shape[0]):
        for j in range(96*VW.shape[1],(96+60)*VW.shape[1]):
            ZZ[i][j-96*VW.shape[1]] = Z[i][j]
    for i in x:
        theta, phi = single_rotation(VW,Z,L,index=i)
        x=np.arange(theta.shape[0])
        pyplot.plot2(x,theta,x,phi,label1='theta',
                label2='phi',save='sphere_rotation%i'%i,showleg=True)
    sum_t = [0.0 for i in range(VW.shape[0]-1)]
    sum_p = [0.0 for i in range(VW.shape[0]-1)]
    for k in range(VW.shape[1]):
        theta, phi = single_rotation(VW,Z,L,index=k)
        for i in range(len(theta)-1):
            x = abs(theta[i+1]-theta[i])
            if x>6/2.0:
                x = abs(x - 6.14)
            sum_t[i] += x
        for i in range(len(phi)-1):
            y = abs(phi[i+1]-phi[i])
            if y>3/2.0:
                y = abs(y - 3.14)
            sum_p[i] += y
    for i in range(len(sum_t)-1):
        sum_t[i+1] += sum_t[i]
        sum_p[i+1] += sum_p[i]
    for i in range(len(sum_t)):
        sum_t[i] = sum_t[i]/VW.shape[1]
        sum_p[i] = sum_p[i]/VW.shape[1]
    x=range(len(sum_p))
    pyplot.plot2(x,sum_t,x,sum_p,label1='theta',
            label2='phi',save='sphere_total_rotation%i'%i,showleg=True)


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
def structure_factor(VW, L, var, n_frames=10,
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
def species_mixing(V,W,L,n_finish=1,n_start=0,delta=100,rcut=20):
    #The total number of frames we are going to look at
    if rcut == False:
        rcut=min_particle_distance(V[-1:],W[-1:],L)+4
    print rcut
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
    pyplot.plot2(x, s, x, o, label1='Same',
            label2='other', save='mixing_avg', showleg=True)
#####################################################
# Look at defects 
#####################################################
def find_defects(V,W,L,n_finish=1,n_start=0,delta=20,rcut=False):
    #The total number of frames we are going to look at
    if rcut == False:
        rcut=min_particle_distance(V[-2:-1],W[-2:-1],L)
    n_frames = (n_finish-n_start)/delta
    x = np.arange(n_start,n_finish,delta)
    defects_V = np.zeros((len(x),15))
    defects_W = np.zeros((len(x),15))
    # Find the number of defects in each frame
    count = 0
    for k in x:
        print k
        for i in range(V[k].shape[0]):
            #find the points that are nearest neighbor that are different
            defects_W[count] += len(nearest_neighbors_point(W[k],W[k][i],L,rcut)[0])
            defects_V[count] += len(nearest_neighbors_point(V[k],V[k][i],L,rcut)[0])
        defects_W[count] = defects_W[count]/8
        defects_V[count] = defects_V[count]/8
        count+=1
    print defects_W
    print defects_V
    pyplot.plot2(x, defects_W, x, defects_V, save='defects_per_frame')
## \brief read mylog file and makes a plot of the energy
#
# \returns x coordinates of time step
# \returns y whatever value of the row you are looking for
#
# \param row the row you would like to grab data from 
def mylog(row = 3):
    import re
    #read log file
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

#################################################
# What to run
#################################################
def run_all():
    #setup
    MD.util.dcd_to_xyz()
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    last = MD.util.get_total_frames('dna')
    var = MD.util.get_path_variables()
    L = MD.util.get_box_size()
    M=MD.ReadCord(frames=last)
    V=M.cord_auto(['V'])
    W=M.cord_auto(['W'])
    VW=M.cord_auto(['V','W'])
    VW,V,W = drift_remove_all(VW,V,W,L)
    rcut=False
    delta=50
    #analysis
    try:
        print 'Finding MSD'
        msd(VW,L)
        plt.close()
    except:
        print 'MSD_end failed'
    try:
        print "finding g(s)"
        distance_distribution(V,W,L,n_start=last-20)
        plt.close()
    except:
        print 'g(s) failed'
    try:
        print 'finding networks'
        find_networks(M, VW, L, var, last, delta=delta)
        plt.close()
    except:
        print 'networks failed'
    try:
        print 'finding solid particles'
        solid_particles_vmd(VW,L,var,rcut=rcut)
        plt.close()
    except:
        print 'solid particles failed'
    try:
        print 'finding cube angle'
        cube_angle_connected(M,VW,L,var['ndna'])
        plt.close()
    except:
        print 'cube angle failed'
    try:
        print 'finding rotations'
        rotations(M,VW,L)
        plt.close()
    except:
        print 'rotations failed'
    try:
        print 'finding sf'
        structure_factor(VW, L, var, n_start=last-20)
        plt.close()
    except:
        print 'sf failed'
    try:
        print 'Finding Speceies Mixing'
        species_mixing(V,W,L,n_finish=last,delta=delta,rcut=rcut)
        plt.close()
    except:
        print 'species mixing failed'
    try:
        find_defects(V,W,L,n_finish=last,delta=delta)
        plt.close()
    except:
        print 'finding defects failed'
    try:
        find_lifetime(M,L,n_finish=last,n_start=last-40)
    except:
        print 'finding lifetime failed'
    try:
        mylog()
    except:
        print 'mylog failed'
def run_debug():
    #setup
    MD.util.dcd_to_xyz()
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    last = MD.util.get_total_frames('dna')
    var = MD.util.get_path_variables()
    L = MD.util.get_box_size()
    M=MD.ReadCord(frames=last)
    V=M.cord_auto(['V'])
    W=M.cord_auto(['W'])
    VW=M.cord_auto(['V','W'])
    VW,V,W = drift_remove_all(VW,V,W,L)
    delta=60
    #analysis
    print 'Finding MSD'
    msd(VW,L)
    plt.close()
    print "finding g(s)"
    distance_distribution(V,W,L,n_start=last-20)
    plt.close()
    print 'finding networks'
    find_networks(M, VW, L, var, n_finish=last,delta=delta)
    plt.close()
    print 'networks failed'
    solid_particles_vmd(VW,V,W,L,var)
    plt.close()
    print 'finding sf'
    structure_factor(VW, L, var, n_start=last-20)
    plt.close()
    print 'Finding Speceies Mixing'
    species_mixing(V,W,L,n_finish=last)
    plt.close()
    print 'finding defects'
    find_defects(V,W,L,n_finish=last)
    print 'find lifetime'
    find_lifetime(M,L,n_finish=last,n_start=last-40)
def run_single():
    #setup
    MD.util.dcd_to_xyz()
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    last = MD.util.get_total_frames('dna')
    var = MD.util.get_path_variables()
    L = MD.util.get_box_size()
    M=MD.ReadCord(frames=last)
    V=M.cord_auto(['V'])
    W=M.cord_auto(['W'])
    VW=M.cord_auto(['V','W'])
    VW,V,W = drift_remove_all(VW,V,W,L)
    #quick fix for a bug
    delta=last/7
    #analysis
    #analysis
    msd(VW,L)
    plt.close()
    structure_factor(VW, L, var, n_start=last-20)
    plt.close()
    distance_distribution(V,W,L,n_start=last-20)
    plt.close()
    #cube_angle_connected(M, L, var['ndna'], frame=-1)
    #plt.close()
    find_networks(M, VW, L, var, last,delta=delta)
    plt.close()
    #find_lifetime(M,L,n_finish=last,n_start=last-40)
    #plt.close()
    solid_particles_vmd(VW,V,W,L,var,skip=delta)
    plt.close()
    #mylog()
    #plt.close()
    species_mixing(V,W,L,n_finish=last)
    plt.close()
    rotations(M,V,L,'Z')
    plt.close()
    #cube_angle(M,L)
    #plt.close()
    #not_connected_nn(M,VW, L, var['ndna'], n_finish=last, n_start=0, delta=delta)
    #plt.close()
    find_defects(V,W,L,n_finish=last,n_start=0,delta=delta)
def run_test():
    #setup
    MD.util.dcd_to_xyz()
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    last = MD.util.get_total_frames('dna')
    var = MD.util.get_path_variables()
    L = MD.util.get_box_size()
    M=MD.ReadCord(frames=last)
    try:
        V = util.pickle_load('V.pkl')
        W = util.pickle_load('W.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
        V=M.cord_auto(['V'])
        W=M.cord_auto(['W'])
        VW=M.cord_auto(['V','W'])
        VW,V,W = drift_remove_all(VW,V,W,L)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(W,'W.pkl')
        util.pickle_dump(VW,'VW.pkl')
    x = range(0,last-20,last/3)
    x.append(last-20)
    print x
    #mylog()
    #plt.close()
    #solid_particles_vmd(VW,V,W,L,var,skip=50)
    #plt.close()
    species_mixing(V,W,L,n_finish=last,rcut=28)
    plt.close()
    #for i in x:
    #    print "finding g(s)"
    #    structure_factor(VW, L, var, n_start=i,save='sf'+str(i))
    #    distance_distribution(V,W,L,n_start=i)
    #    plt.close()
    #msd(VW,L)
    #plt.close()
    #end_end(V,W,M,L,var)
    #find_networks(M, VW, L, var, last,delta=50)
    #plt.close()
    #quick fix for a bug
    #rotations(M,V,L,'Z')
    #sphere_rotations(M,W,L,'N')
    #rotations(M,VW,L,'Z')
    #cube_second_angle(M,L)
    #cube_connected_face(M, L, var['ndna'], VW.shape[1], n_finish = last, delta = 30)
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

