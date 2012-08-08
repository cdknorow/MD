# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MD
############ 
# Run from inside a folder
#########
reload(MD)
reload(pyplot)
#################################################
# Number of Connections
#################################################
def find_clusters(M,L,var,n_finish,n_start=1,delta=30,rcut=1.0,step=5e4):
    B=M.cord_auto(['B'])
    #The total number of frames we are going to look at
    n_frames = (n_finish-n_start)/delta
    x=np.arange(n_start,n_finish,delta)
    print x
    steps=x*step
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    try:
        connections=util.pickle_load('con.pkl')
    except:
        connections = con.connections(C[x],G[x],L,rcut=rcut)
        util.pickle_dump(connections,'con.pkl')
    connections=con.num_connections(connections,var['nsphere'])
    pyplot.plot(steps,connections,xlabel='Time',
            ylabel='hyrbid. density', save='connections')
    return steps,connections
#####################################################
#find number of a specified pair of particles that are within r_cut of each
#####################################################
def find_networks(M, L, var, n_finish, n_start=1, delta=30, rcut=1.0):
    C=M.cord_auto(['C'])
    G=M.cord_auto(['G'])
    #The total number of frames we are going to look at
    x=np.arange(n_start,n_finish,delta)
    #Find the number of connections at specific points x
    import MD.analysis.connections as con
    import MD.analysis.percolate as perc
    try:
        connections = util.pickle_load('con.pkl')
    except:
        connections = con.connections(C[x],G[x],L,rcut=rcut)
        util.pickle_dump(connections,'con.pkl')
        connections = con.num_connections(connectsion,var['Nsphere'])
    percolated, networks = perc.percolate(connections,var['ndna'],var['nsphere'])
    pyplot.plot(x,networks,xlabel='Time',
            ylabel='networks', save='networks')
    return x, networks, connections
#####################################################
# Find Static Structure Facto
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor(VW, L, var, n_frames=10, n_start=1, 
        filters=0.05, dh=0.05):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    stmp,qx = sf.sfactor(VW[n_start:n_start+n_frames],L=L,l=10,
            N_points=var['nsphere'],n_frames=n_frames) 
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx,dh=dh)
    #Plot the graph
    xlabel = '$|\\vec{q}$ $|$'
    ylabel = '$S(\\vec{q}$ $)$ '
    save = 'sf'
    pyplot.plot(Q, S, linestyle='', marker='x', 
            xlabel=xlabel, ylabel=ylabel, save=save)
    QQ,SS,PM = sf.sf_filter(Q,S,primitive_vectors,filters=0.05)
    if len(QQ)>500:
        bars,height = sf.sf_max(QQ,qpeaks=8)
        db_x,db_y = sf.dbfactor(msdr=5)
        pyplot.plot_bar(bars, height, xlabel=xlabel, ylabel=ylabel,
                save=save)
        pyplot.plot(db_x, db_y, linestyle='--', xlabel=xlabel, ylabel=ylabel,
                save=save)
        return Q, S, bars, height, db_x, db_y
    return Q, S, primitive_vectors

#################################################
# What to run
#################################################
def run_all():
    #setup
    MD.util.dcd_to_xyz()
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    last = MD.util.get_total_frames('nrod')
    var = MD.util.get_path_variables()
    L = MD.util.get_box_size()
    M=MD.ReadCord(frames=last)
    B=M.cord_auto(['B'])
    try:
        print "Finding Miscell Clusters"
        find_cluster(B,L)
    try:
        print 'Finding MSD'
        msd(B,L)
        plt.close()
    except:
        print 'MSD_end failed'
    try:
        print 'finding connections'
        find_connections(M,L,var,n_finish=last,delta=30)
        plt.close()
    except:
        print 'connections failed'
    try:
        print 'finding networks'
        find_networks(M, L, var, last)
        plt.close()
    except:
        print 'networks failed'
    try:
        print 'finding sf'
        structure_factor(VW, L, var, n_start=last-20)
    except:
        print 'sf failed'
#For multiple directories
if __name__ == '__main__':
    for f in sorted(os.listdir("./")):
        if os.path.isdir(f):
            os.chdir(f)
            try:
                run_all()
            except:
                print "Someting Failed"
            print '##\nfinished with:',f
            os.chdir('../')
#For single directories
#if __name__ == '__main__':
#    run_all()
