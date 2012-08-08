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
############################################
## returns the arrays with the drift removed
############################################
def drift_remove(VW,L,index):
    from MD.analysis.drift_remove import eliminate_drift
    VW = eliminate_drift(VW,L,index)
    return VW
    #Find the msd of the system
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
    sp.crystal_size(L,VW.shape[0]/skip,VW.shape[1])
######################
# FInd the Angle between cube normal and c-axis
# caxis is the vector of the caxis
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
def run():
    #setup
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    print L
    last = M.frames
    try:
        V = util.pickle_load('V.pkl')
    except:
        V=M.cord_auto(['V'])
        V_index=M.get_names(['V'])
        V = drift_remove(V,L,V_index)
        util.pickle_dump(V,'V.pkl')
    x = range(0,last,last/8)
    #mylog()
    #plt.close()
    #for i in x:
    #    #print "finding s(q)"
    #    structure_factor(V, L, var, n_start=i,save='sf'+str(i))
        #print "finding g(s)"
        #distance_distribution(V,V,L,n_start=i)
        #plt.close()
    msd(V,L)
    plt.close()
    #find_lifetime(M,L,n_finish=last,n_start=last-30,rcut=1.0)
    #plt.close()
    #find_networks(M, V, L, var, last,delta=10)
    #plt.close()
    #plt.close()
    #quick fix for a bug
#For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            try:
#                run()
#                print '##\nfinished with:',f
#            except:
#                print 'failed in this directory'
#            os.chdir('../')
#For single directories
if __name__ == '__main__':
    run()

