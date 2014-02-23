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
#get connections
def _connections(M):
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
def find_networks(M, VW, L,x,ndna=25):
    import MD.canalysis.connections as con
    import MD.analysis.connections as conn
    import MD.analysis.graph as graph
    #plot total number of connections
    connections = _connections(M)
    num_connections=conn.num_connections(connections,VW.shape[1])
    con_all = []
    for i in range(len(connections)):
        con_all.append(len(connections[i][0]))
    print len(x),len(num_connections)
    pyplot.plot(x,num_connections,xlabel='Time',
            ylabel='hyrbid. density', save='connections')
    pyplot.plot(x,con_all,xlabel='Time',
            ylabel='hyrbid', save='connections_all')
    plt.close()
    #get the info
    networks, num_networks, deg, neighbors, num_n, gr = graph.grapher(connections,VW.shape[1],ndna)
    util.pickle_dump(networks,'net.pkl')
    #plot the number of neighbors at each timesteps
    pyplot.plot(x,num_n,xlabel='Time',
            ylabel='num neighbors', save='neighbors')
    plt.close()
    print 'making plot'
    net = []
    for i in networks:
        net.append(len(i))
    pyplot.plot(x,net,xlabel='t',ylabel='networks',save='net')
    label = ['networks','1','2','3','4','5','6','7','8','9']
    x = [x for i in range(len(deg))]
    #pyplot.plot(x,num_networks,xlabel='t',ylabel='Networks',save='net')
    #pyplot.plot_multi(x[:4],deg[:4],label,xlabel='time',ylabel='number',save='con_net')
    plt.close()
    return networks, num_networks, deg, neighbors, num_n, gr

    #pyplot.plot_multi(x[:4],deg[:4],label,xlabel='time',ylabel='number',save='con_net')
    plt.close()
    return networks, num_networks, deg, neighbors, num_n, gr
def find_lifetime(M):
    import MD.analysis.lifetime as lifetime
    conn = _connections(M)
    life = []
    for i in range(3):
        life.append(lifetime.lifetime(conn[i:])[:len(conn)-3])
    life_avg = []
    for i in range(len(life[0])):
        avg = 0
        for j in range(len(life)):
            avg += life[j][i]/float(len(life))
        life_avg.append(avg)
    pyplot.plot(range(len(life_avg)),life_avg,xlabel='Time',
            ylabel='remain connections', save='lifetime')
def center_end_connected(M,L):
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
    con = _connections(M)
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
def end_end_connected(M,L): #find the length of the polymer
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
    con = _connections(M)
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
def polymer_gauss(M,VW, Z,L, frames, rcut=1.0, step=5e4):
    #\brief find the gauss map of polymers surrounding NC
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
    util.print_box_volume(L_cont,delta=delta)
    last = M.frames
    #V_index = M.get_index(['V'])
    #W_index = M.get_index(['W'])
    try:
        VW = util.pickle_load('VW.pkl')
    except:
         VW=M.cord_auto(['V','W'])
         #VW,Z,W = drift_remove_all(VW,Z,VW,L)
         util.pickle_dump(VW,'VW.pkl')
    #L_last = L_cont[0][0]
    #for i,j in enumerate(L_cont):
    #    if j[0] != L_last:
    #        L_last = j[0]
    #        if i - 5 > 0:
    #            x.append(i-5)
    x = range(0,last)
    index = M.get_index(['C','G'])
    ndna = len(index)/VW.shape[1]
    #find_networks(M, VW, L_cont,x,ndna= ndna)
    find_lifetime(M)
#################################################
# dcd
#################################################
def ipy_single():
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
    util.print_box_volume(L_cont,delta=delta)
    last = M.frames
    #V_index = M.get_index(['V'])
    #W_index = M.get_index(['W'])
    try:
        VW = util.pickle_load('VW.pkl')
    except:
         VW=M.cord_auto(['V','W'])
         #VW,Z,W = drift_remove_all(VW,Z,VW,L)
         util.pickle_dump(VW,'VW.pkl')
    #L_last = L_cont[0][0]
    #for i,j in enumerate(L_cont):
    #    if j[0] != L_last:
    #        L_last = j[0]
    #        if i - 5 > 0:
    #            x.append(i-5)
    x = range(0,last)
    index = M.get_index(['C','G'])
    ndna = len(index)/VW.shape[1]
    return [M,VW,L_cont,x,ndna]
##For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            run_single()
#            print '##\nfinished with:',f
#            os.chdir('../')
#For single directories
if __name__ == '__main__':
    run_single()

