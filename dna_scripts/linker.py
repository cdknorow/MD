# -*- coding: utf-8 -*-
import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
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
#get connections
def connections(M):
    try:
        con = util.pickle_load('conn.pkl')
    except:
        fid = open('conn_who.dat','r')
        c_index = M.get_index(['C'])
        g_index = M.get_index(['G'])
        cc_index = M.get_index(['CC'])
        gg_index = M.get_index(['GG'])
        con = []
        count = 0
        for line in fid.readlines():
            con_K = []
            con_T = []
            print line
            #linker alternates between c g and cc gg
            if count%2 == 0:
                line = line.replace('}','')
                for i in line.split('{')[1].split():
                    print i
                    #translate connections
                    con_K.append(c_index.index(int(i)))
                for i in line.split('{')[2].split():
                    #translate connections
                    con_T.append(g_index.index(int(i)))
            else:
                line = line.replace('}','')
                for i in line.split('{')[1].split():
                    #translate connections
                    con_K.append(cc_index.index(int(i)))
                for i in line.split('{')[2].split():
                    #translate connections
                    con_T.append(gg_index.index(int(i)))
            con.append([con_K,con_T])
            count += 1
        del con[-1]
        util.pickle_dump(con,'conn.pkl')
        print con
        return con
#\brief find the gauss map of connections surrounding NC
def linker_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    A = 1
    ndna=150
    C = M.cord_auto(['C'])
    CC = M.cord_auto(['CC'])
    try:
        conn = util.pickle_load('conn.pkl')
    except:
        conn = connections(M)
        util.pickle_dump(conn,'conn.pkl')

    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    for k in frames:
        print k
        conn_location = []
        for i in conn[k*2][0]:
            #we must rotate about a specific cube reference frame
            cube = i/ndna
            V = VW[k][cube]
            x_r = points.unit(points.vector1d(V,Z[k][cube*6],L))
            y_r = points.unit(points.vector1d(V,Z[k][cube*6+2],L))
            z_r = points.unit(points.vector1d(V,Z[k][cube*6+4],L))
            v = np.array([x_r,y_r,z_r])
            R = points.reference_rotation(v,w)
            d = points.dist(V,C[k][i],L)[0]
            c_r = points.unit(points.vector1d(V,C[k][i],L))
            conn_location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
        for i in conn[k*2+1][0]:
            #A is the number of A cubes
            cube = i/ndna
            V = VW[k][cube]
            x_r = points.unit(points.vector1d(V,Z[k][cube*6],L))
            y_r = points.unit(points.vector1d(V,Z[k][cube*6+2],L))
            z_r = points.unit(points.vector1d(V,Z[k][cube*6+4],L))
            v = np.array([x_r,y_r,z_r])
            R = points.reference_rotation(v,w)
            c_r = points.unit(points.vector1d(V,CC[k][i],L))
            conn_location.append(points.unit(np.dot(R,np.transpose(c_r))))
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
def test_linker_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    A = 1
    ndna=100
    C = M.cord_auto(['C'])
    CC = M.cord_auto(['CC'])

    gauss_map = []
    w = np.array([[1,0,0],[0,1,0],[0,0,1]])
    test_k = np.array([[0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],
                        [-0.5,0.5,0.5],[-0.5,-0.5,0.5]])
    for k in frames:
        conn_location = []
        for i in test_k:
            #we must rotate about a specific cube reference frame
            cube = 0
            V = VW[k][cube]
            x_r = points.unit(points.vector1d(V,Z[k][cube*6],L))
            y_r = points.unit(points.vector1d(V,Z[k][cube*6+1],L))
            z_r = points.unit(points.vector1d(V,Z[k][cube*6+5],L))
            v = np.array([x_r,y_r,z_r])
            R = points.reference_rotation(v,w)
            print R
            print v
            print '######'
            print np.dot(R,v)
            print points.unit(np.dot(R,i))
            conn_location.append(points.unit(np.dot(R,i)))
        gauss_map.append(conn_location)
    fid = open('gaussmap_connections_test.xyz','w')
    max_gauss = 0
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss,frames[k]))
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in range(max_gauss - len(gauss_map[k])):
            fid.write('N  0 0 0\n')
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
############################################
# lifetime of connections
########################################
def find_lifetime(M,L,steps,step_range=30,delta=4,rcut=1.0,step=5e4):
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
        print 'Frame',i
        x=np.arange(i,step_range+i,delta)
        #Find he number of connections at specific points x
        remain = life.lifetime(C[x],G[x],L)
        x=np.arange(i,step_range+i,delta)*50000
        print remain
        pyplot.plot(x,remain,xlabel='time', ylabel='remaining connections',
                save='lifetime%i'%i)
        plt.close()
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
def print_box_volume(L_cont,delta=10):
    fid = open('box_length.txt','w')
    T = 5+8**0.6
    fid.write('frame  L   phi\n') 
    for i in range(0,len(L_cont),delta):
        if L_cont[i][0] != 0:
            phi = 54*T**3/L_cont[i][0]**3
            fid.write(('%i   %.2f   %.2f\n')%(i,L_cont[i][0],phi))
    fid.close()
#################################################
# What to run
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
    #drift_remove_all(VW,Z,VW,L)
    #x = range(0,last,delta)
    x = [100]
    #cubic_order_animate(VW,Z,L_cont,x)
    x = range(400,last)
    linker_gauss(M,VW,Z, L, x)
    #test_linker_gauss(M,VW,Z, L, x)
    #cubic_order(VW,Z,L_cont,x)
    #cube_connected_rotate(VW,Z,M, L_cont)
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
    #rotation_diffusion_split(Z,VW,L_cont,x)
    import MD.analysis.cubatic_order as c
    print VW.shape
    #c.gaussmap(VW,Z,L_cont,x,3)
#################################################
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
    run_single()
