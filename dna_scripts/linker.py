# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
from MD.analysis.particle_distance import particle_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.plot.histogram import histogram_normal as histogram
from MD.analysis.rotation import single_rotation
############ 
# Run from inside a folder
#########
reload(MD)
reload(pyplot)
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
    pyplot.plot(x, y, xlabel='time', ylabel=label, save=save)
    return x, y
#####################################################
# Look at mixing of species of spheres
#####################################################
def species_mixing_initial(V,W,L,frame,rcut=False):
    from MD.analysis.particle_distance import min_particle_distance
    if rcut == False:
        rcut=min_particle_distance(V[-1:],W[-1:],L)+5
    k=frame
    Same = np.zeros((V[0].shape[0],1))
    Other = np.zeros((V[0].shape[0],1))
    for i in range(V[frame].shape[0]):
        Same[i] = len(nearest_neighbors_point(V[k],V[k][i],L,rcut=rcut)[0])
        Other[i] = len(nearest_neighbors_point(W[k],V[k][i],L,rcut=rcut)[0])
    x = np.arange(0, Same.shape[0])
    pyplot.plot2(x, Same, x, Other, label1='Same', label2='other',
            save='mixing_single', showleg=True)
#####################################################
# Look at mixing of species of spheres
#####################################################
def species_mixing(V, W, L, n_finish = 1, n_start = 0, delta = 20, rcut = False):
    #The total number of frames we are going to look at
    from MD.analysis.particle_distance import min_particle_distance
    if rcut == False:
        rcut=min_particle_distance(V[-1:],W[-1:],L)+1
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
    pyplot.plot2(x, s, x, o, label1='Same',
            label2='other', save='mixing_avg', showleg=True)
#####################################################
# Look at defects 
#####################################################
def find_defects(V,W,L,n_finish=1,n_start=0,delta=20,rcut=False):
    #The total number of frames we are going to look at
    from MD.analysis.particle_distance import min_particle_distance
    if rcut == False:
        rcut=min_particle_distance(V[-1:],W[-1:],L)+1
    n_frames = (n_finish-n_start)/delta
    x = np.arange(n_start,n_finish,delta)
    defects = np.zeros((len(x),15))
    # Find the number of defects in each frame
    count = 0
    for k in x:
        print k
        for i in range(V[k].shape[0]):
            #find the points that are nearest neighbor that are different
            defects[count] += len(nearest_neighbors_point(W[k],W[k][i],L,rcut)[0])
            defects[count] += len(nearest_neighbors_point(V[k],V[k][i],L,rcut)[0])
        count+=1
    pyplot.plot(x, defects, save='defects_per_frame')
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
    stmp,qx = sf.sfactor(VW[n_start:n_start+n_frames],L=L,l=10)
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
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
############################################
## 	g(r)
############################################
def distance_distribution(V,W,L,n_frames=15,n_start=1):
    #finds end to end distances and plots a histogram with that data
    def hist(M1,M2,save_name,label):
        print M1
        print M2
        distance=particle_distance(M1,M2,L)
        print distance
        hist_s,xs,max_hist=histogram(distance,bins=50)
        return xs,hist_s
    ########################################
    #Plot Data for AA,AB,BB distances
    save_name='_time_'+'%i'%(n_start)
    start=n_start
    finish=n_start+n_frames
    max_hist=0
    AB=hist(V[start:finish],W[start:finish],save_name,'A-B')
    AA=hist(W[start:finish],W[start:finish],save_name,'A-A')
    BB=hist(V[start:finish],V[start:finish],save_name,'B-B')
    pyplot.plot3(AA[0],AA[1],AB[0],AB[1],BB[0],BB[1],'s','g(s)',
            'A-A','A-B','B-B',save='nnplot_'+save_name,showleg=True)
############################################
## msd
############################################
def msd(VW,L,step=5e4, save = 'MSD_time'):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x,msd=msd(VW,L,step=step)
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save=save)
######################
# Find the solid particles in a simulation
#####################
def solid_particles(VW, L, var, skip=30, step=5e4, rcut = False):
    from MD.analysis import bond_order
    from MD.analysis.particle_distance import min_particle_distance
    if rcut == False:
        rcut=min_particle_distance(VW[-1:],VW[-1:],L)+6
    BCC = []
    SC = []
    try:
        BCC = util.pickle_load('crystalbcc.pkl')
        SC = util.pickle_load('crystalsc.pkl')
    except:
        for i in range(0,len(VW),skip):
            print 'step', i
            BCC.append(bond_order.structure_correlation(VW[i], L, c_cut=0.1, rcut=rcut, l=6))
            SC.append(bond_order.structure_correlation(VW[i], L, rcut = rcut, l=4,
                c_cut = 0.1, crystal = 4))
            util.pickle_dump(BCC, 'crystalbcc.pkl')
            util.pickle_dump(SC, 'crystalsc.pkl')

    x = np.arange(len(BCC))
    x *= skip * step
    xlabel = 'time'
    ylabel = 'f(solid)'
    SC = np.array(SC) / float(var['nsphere'])
    BCC = np.array(BCC) / float(var['nsphere'])
    pyplot.plot2(x, SC, x, BCC, label1 = 'SC', label2 = 'BCC', xlabel = xlabel,
            ylabel = ylabel, save = 'crystals', showleg = True)

    return x, BCC, SC
#########################################
# What you want to do in those directories
########################################
def end_end(M,L,var):
    #find the length of the polymer
    KT=M.cord_auto(['K','TT'])
    S=M.cord_auto(['M'])
    edge=[]
    print L
    for k in range(S.shape[0]):
        for i in range(S.shape[1]):
            edge.append(points.dist(KT[k][i],S[k][i],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,50)

    pyplot.plot_bar(xedge,hist_edge,save='edge_distance')
    plt.close()
    #find the length of the linker
    TF=M.cord_auto(['G'])
    TL=M.cord_auto(['CC'])
    edge=[]
    for k in range(TF.shape[0])[-1:]:
        for i in range(TF.shape[1]):
            edge.append(points.dist(TF[k][i],TL[k][i],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,50)

    pyplot.plot_bar(xedge,hist_edge,save='linker_distance')
    plt.close()
## \brief Number of Connections betweeen linkers
#
# \returns steps x frames number for each data point
# \returns connections_A,connections_B list of connections in each frame
# \returns one, two number of connections in each frame
#
# \param M lets you read variables into a Matrix ie M.courd_auto(['G'])
# \param L Length of box
# \param var dictionary of variables read from the file path name
# \param n_finish the last frame to look at
# \param n_start is the first frame to look at
# \param delta the number of frames to skip over
# \param rcut the distance a connection is defined by 
# \param step lets the plot know have many timesteps between each frame
#
def find_linker_connections(M, L, var,n_finish, n_start=1,
        delta=30, rcut=1.0, step=5e4):
    #The total number of frames we are going to look at
    n_frames = (n_finish - n_start) / delta
    x = np.arange(n_start,n_finish,delta)
    print x
    steps = x * step
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    reload(con)
    try:
        connections_A = util.pickle_load('conA.pkl')
        connections_B = util.pickle_load('conB.pkl')
    except:
        try:
            C = util.pickle_load('C.pkl')
            G = util.pickle_load('G.pkl')
            CC = util.pickle_load('CC.pkl')
            GG = util.pickle_load('GG.pkl')
        except:
            C = M.cord_auto(['C'])
            G = M.cord_auto(['G'])
            CC = M.cord_auto(['CC'])
            GG = M.cord_auto(['GG'])
            C = util.pickle_dump(C,'C.pkl')
            G = util.pickle_dump(G,'G.pkl')
            CC = util.pickle_dump(CC,'CC.pkl')
            GG = util.pickle_dump(GG,'GG.pkl')
        #put the linker second to keep track
        connections_A = con.connections(C[x],G[x],L,rcut=rcut)
        connections_B = con.connections(GG[x],CC[x],L,rcut=rcut)
        util.pickle_dump(connections_A,'conA.pkl')
        util.pickle_dump(connections_B,'conB.pkl')
    one,two,l1,l2 = con.both_connections(connections_A,connections_B)
    connections_A = con.num_connections(connections_A,var['nsphere'])/2.
    connections_B = con.num_connections(connections_B,var['nsphere'])/2.
    pyplot.plot2(steps, connections_A, steps, connections_B,
            label1='C',label2='CC', xlabel='time', ylabel='hyrbid. density',
            save='connections',showleg=True)


    pyplot.plot2(steps, np.array(one)/27., steps, np.array(two)/27.,
            label1='One Connection',label2='Two Connections',
            xlabel='time', ylabel='hybrid density', save='link_con_hy',
            showleg=True)
    util.pickle_dump([one,two],'one_two.pkl')
    #three = [CC.shape[1] for i in range(len(one))]
    three = [800 for i in range(len(one))]
    #print 'linkers'
    #print CC.shape[1]
    three = np.subtract(three, one)
    three = np.subtract(three, two)
    pyplot.plot3(steps, one, steps, two, steps, three,
            label1='One Connection',label2='Two Connections', label3='No Connections',
            xlabel='time', ylabel='num con', save='link_con',
            showleg=True)
    plt.close()
    return steps, connections_A, connections_B, one, two
def find_linker_fcc_connections(M, L, var,n_finish, n_start=1,
        delta=30, rcut=1.0, step=5e4):
    #The total number of frames we are going to look at
    n_frames = (n_finish - n_start) / delta
    x = np.arange(n_start,n_finish,delta)
    print x
    steps = x * step
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    reload(con)
    try:
        connections_A = util.pickle_load('conA.pkl')
        connections_B = util.pickle_load('conB.pkl')
    except:
        C = M.cord_auto(['C'])
        G = M.cord_auto(['G'])
        #put the linker second to keep track
        connections_A = con.connections(C[x],G[x][range(1,G.shape[1],2)],L,rcut=rcut)
        connections_B = con.connections(C[x],G[x][range(0,G.shape[1],2)],L,rcut=rcut)
        util.pickle_dump(connections_A,'conA.pkl')
        util.pickle_dump(connections_B,'conB.pkl')
    one,two,l1,l2 = con.both_connections(connections_A,connections_B)
    connections_A = con.num_connections(connections_A,var['nsphere'])
    connections_B = con.num_connections(connections_B,var['nsphere'])
    pyplot.plot2(steps, connections_A, steps, connections_B,
            label1='C',label2='CC', xlabel='time', ylabel='hyrbid. density',
            save='connections')

    #three = [CC.shape[1] for i in range(len(one))]
    three = [800 for i in range(len(one))]
    #print 'linkers'
    #print CC.shape[1]
    three = np.subtract(three, one)
    three = np.subtract(three, two)
    pyplot.plot3(steps, one, steps, two, steps, three,
            label1='One Connection',label2='Two Connections', label3='No Connections',
            xlabel='time', ylabel='num con', save='link_con',
            showleg=True)
    plt.close()
    return steps, connections_A, connections_B, one, two
def find_linker_distance(M, L, var, frame = -1,rcut=1.0, step=5e4):
    #The total number of frames we are going to look at
    #Find he number of connections at specific points x
    import MD.analysis.connections as con
    reload(con)
    C = M.cord_auto(['C'])
    G = M.cord_auto(['G'])
    CC = M.cord_auto(['CC'])
    GG = M.cord_auto(['GG'])
    #put the linker second to keep track
    connections_A = con.connections(C[frame:],G[frame:],L,rcut=rcut)
    connections_B = con.connections(GG[frame:],CC[frame:],L,rcut=rcut)
    print connections_A
    print connections_B
    one,two,l1,l2 = con.both_connections(connections_A,connections_B)
    #Get the e2e distance for the ones
    edge=[]
    bins=13
    for i in l1[-1]:
        edge.append(points.dist(G[frame][i-C.shape[1]],CC[frame][i-C.shape[1]],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,bins)
    plt.close()
    pyplot.plot_bar(xedge,hist_edge,save='linker_1con_distance')
    plt.close()
    edge=[]
    for i in l2[-1]:
        edge.append(points.dist(G[frame][i-C.shape[1]],CC[frame][i-C.shape[1]],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,bins)
    pyplot.plot_bar(xedge,hist_edge,save='linker_2con_distance')
    plt.close()
    #get a list of linkers without any connections
    l0 = range(C.shape[1],C.shape[1]+CC.shape[1])
    l0 = points.difference(l0, l1[-1])
    l0 = points.difference(l0, l2[-1])
    print len(l0)
    print len(l1[0])+len(l2[0])+len(l0)
    edge=[]
    for i in l0:
        edge.append(points.dist(G[frame][i-C.shape[1]],CC[frame][i-C.shape[1]],L)[0])
    hist_edge,xedge,max_edge=histogram(edge,bins)
    pyplot.plot_bar(xedge,hist_edge,save='linker_0con_distance')
    plt.close()
    from MD.analysis.msd import msd
    x,msd=msd(CC,L,step=50000)
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save='msd_linker')
    return edge
## \brief find number of a specified pair of particles that are within r_cut of each
#
# \returns x cordinates of x graph
# \returns networks a graph containing the independent sets in each frame
# \returns num_networks array containing the number of independent sets in each frame
#
# \param M lets you read variables into a Matrix ie M.courd_auto(['G'])
# \param VW matrix of V and W elements
# \param L Length of box
# \param var dictionary of variables read from the file path name
# \param n_finish the last frame to look at
# \param n_start is the first frame to look at
# \param delta the number of frames to skip over
# \param rcut the distance a connection is defined by 
#
def find_networks(M, VW, L, var, n_finish, n_start=0, delta=30, rcut=1.0):
    ndna = var['ndna']
    #The total number of frames we are going to look at
    x=np.arange(n_start,n_finish,delta)
    #Find the number of connections at specific points x
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    reload(graph)
    try:
        connections_A = util.pickle_load('NWconA.pkl')
        connections_B = util.pickle_load('NWconB.pkl')
    except:
        C = M.cord_auto(['C'])
        G = M.cord_auto(['G'])
        CC = M.cord_auto(['CC'])
        GG = M.cord_auto(['GG'])
        #put the linker second to keep track
        connections_A = con.connections(C[x],G[x],L,rcut=rcut)
        connections_B = con.connections(GG[x],CC[x],L,rcut=rcut)
        util.pickle_dump(connections_A,'NWconA.pkl')
        util.pickle_dump(connections_B,'NWconB.pkl')
    networks, num_networks = graph.grapher_2(connections_A,
            connections_B, N_nodes = VW.shape[1], N_linker = var['Nlinker'], vertex_points = ndna)
    pyplot.plot(x,num_networks,xlabel='time',
            ylabel='networks', save='networks')
    return x, networks, num_networks
## \brief find number lifetime of connections at a specific frame
#
# \returns x cordinates of x graph
# \returns remains, remains 2 array of number of remaining values in each frame
# afte initial frame
#
# \param M lets you read variables into a Matrix ie M.courd_auto(['G'])
# \param VW matrix of V and W elements
# \param L Length of box
# \param var dictionary of variables read from the file path name
# \param n_finish the last frame to look at
# \param n_start is the first frame to look at
# \param delta the number of frames to skip over
# \param rcut the distance a connection is defined by 
def find_lifetime(M,L,n_finish=30,n_start=1,delta=4,rcut=1.0,step=5e4):
    C=M.cord_auto(['C'])
    G=M.cord_auto(['G'])
    CC=M.cord_auto(['CC'])
    GG=M.cord_auto(['GG'])
    #The total number of frames we are going to look at
    n_frames = (n_finish-n_start)/delta
    x=np.arange(n_start,n_finish,delta)
    print n_start
    print n_finish
    print x
    print C.shape
    print G.shape
    print CC.shape
    print GG.shape
    steps=x*step
    #Find he number of connections at specific points x
    import MD.analysis.lifetime as life
    reload(life)
    remain = life.lifetime(C[x],G[x],L)
    remain2 = life.lifetime(GG[x],CC[x],L)
    print remain
    print remain2
    pyplot.plot2(steps,remain,steps,remain2,xlabel='time',
            ylabel='remaining connections', label1='A Conn', label2='Bconn',save='lifetime')
    return steps,remain,remain2
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
    #analysis
    #try:
    #    print 'Finding End to End'
    #    end_end(V,W,M,L,var)
    #    plt.close()
    #except:
    #    print 'end_end failed'
    try:
        print 'Finding MSD'
        msd(VW,L)
        plt.close()
    except:
        print 'MSD_end failed'
    try:
        print 'finding sf'
        structure_factor(VW, L, var, n_start=last-20)
    except:
        print 'sf failed'
    try:
        print "finding g(s)"
        distance_distribution(V,W,L,n_start=last-20)
        plt.close()
    except:
        print 'g(s) failed'
    try:
        print 'finding connections'
        find_linker_connections(M,L,var,n_finish=last,delta=30)
        plt.close()
    except:
        print 'connections failed'
    try:
        print 'finding newtorks'
        find_networks(M, VW, L, var, last)
    except:
        print 'finding networks failed'
    #try:
    #    print 'finding solid particles'
    #    solid_particles(VW,L,var)
    #    plt.close()
    #except:
    #    print 'solid particles failed'
    #try:
    #    print 'Finding Speceies Mixing'
    #    species_mixing(V,W,L,frame=1)
    #    plt.close()
    #except:
    #    print 'species mixing failed'
    #################################################
def run_quick():
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
    VW = V
    #VW=M.cord_auto(['V','W'])
    #print 'Finding MSD'
    #msd(VW,L)
    #plt.close()
    #print 'finding sf'
    #structure_factor(VW, L, var, n_start=last-20)
    #plt.close()
    ##print "finding g(s)"
    #distance_distribution(V,V,L,n_start=last-20)
    #plt.close()
    #end_end(M,L,var)
    #plt.close()
    find_linker_connections(M,L,var,delta=10,n_finish=last)
    plt.close()
    #find_linker_distance(M, L, var)
#################################################
# What run without exceptions
#################################################
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
    #code to ,,teset
    print 'end2end'
    ##end_end(M,L,var)
    #print 'msd'
    #msd(VW,L)
    #print 'sf'
    #structure_factor(VW, L, var, n_start=last-20)
    #print 'gs'
    #distance_distribution(V,W,L,n_start=last-20)
    print 'con'
    find_linker_connections(M,L,var,n_finish=last,delta=30)
    #print 'net'
    #find_networks(M, VW, L, var, last)
    #print 'sp'
    #solid_particles(VW,L,var)
    #print 'sm'
    #species_mixing(V,W,L,n_finish = last)
    #print 'fd'
    #find_defects(V,W,L,n_finish = last)
    #print 'mylog'
    #mylog(row=3)
    #print 'mylog'
    #mylog(row=2)
    #print 'sp_initial'
    #species_mixing_initial(V,W,L,frame=1)
    #print 'fl'
    #find_lifetime(M,L,n_finish=30,n_start=last-45)
##For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            run_quick()
#            #try:
#            #    run_quick()
#            #except:
#            #    print "Someting Failed"
#            #print '##\nfinished with:',f
#            os.chdir('../')
#For single directories
if __name__ == '__main__':
    run_quick()
