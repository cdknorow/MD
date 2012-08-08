import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
import MD.plot.pyplot_eps as pyplot_eps
import MD.dna_scripts.solidparticle as sp
import MD.dna_scripts.diffusion as diffuse
from MD.analysis.particle_distance import particle_distance
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import close_neighbors_point
from MD.plot.histogram import histogram
from MD.plot.histogram import histogram_reg
from MD.plot.histogram import histogram_normal
from MD.analysis.rotation import single_rotation
from MD.util import readxyz
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
############################################
## returns the arrays with the drift removed
############################################
def drift_remove_all(VW,V,W,L,VW_index):
    from MD.analysis.drift_remove import eliminate_drift
    W = eliminate_drift(W,L)
    V = eliminate_drift(V,L)
    VW = eliminate_drift(VW,L,VW_index)
    return VW,V,W
    #Find the msd of the system
#########################################
# What you want to do in those directories
########################################
def end_end(V,W,M,L):
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
                VK.append(points.dist(V[k][i],K[k][j+i*ndna],L)[0]-3)
                WT.append(points.dist(W[k][i],T[k][j+i*ndna],L)[0]-3)
    hist_K,xK,max_k=histogram(VK,10)
    hist_T,xT,max_T=histogram(WT,10)

    plt.close()
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
############################################ ## msd
############################################
def msd(VW,L,step=1,save='MSD_time',time_scale=1):
    from MD.analysis.msd import msd_no_drift #Find the msd of the system x,msd=msd_no_drift(VW,L,step=step)
    D=0.1
    D2=0.03
    tau=10
    try:
        M = util.pickle_load('msd.pkl')
        x = M[0][0]
        msd = M[1][0]
    except:
        x, msd = msd_no_drift(VW,L)
        x = x*time_scale
        util.pickle_dump([[x],[msd]],'msd.pkl')
    #msd_fit
    fit = []
    for i in x:
        fit.append((2*3*D2*(i) - tau*(1-math.e**(-i/tau)))**0.5)
    msd2 = (6*D*x)**0.5
    #pyplot.plot(x,msd,xlabel='Time',ylabel='msd',save=save)
    #pyplot.plot3(x,msd,x,(6*D*x)**0.5,x,(6*D2*x)**0.5,xlabel='Time',ylabel='msd',label1='msd',label2='%.2f'%D,label3='%.2f'%D2,save=save,showleg=True)
    pyplot.plot3(x,msd,x,(6*D2*x)**0.5,x,msd2,xlabel='Time',ylabel='msd',label1='msd',label2='%.2f'%D2,label3='%.2f'%D,save=save,showleg=True)
def msd_jump_average(VW,L):
    reload(diffuse)
    try:
        jump_d = util.pickle_load('jump_d.pkl')
        jump_t = util.pickle_load('jump_t.pkl')
    except:
        fid = open('large.xyz','r')
        M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
        crystal = M.cord_auto(['V','W'])
        bcc = np.zeros((777,432,1))
        for frame in range(bcc.shape[0]):
            for i in range(bcc.shape[1]):
                if crystal[frame][i][0] > L[0]:
                    bcc[frame][i]=0
                else:
                    bcc[frame][i]=1
        jump_d, jump_t = diffuse.crystal_jump(VW,bcc,L)
        util.pickle_dump(jump_d,'jump_d.pkl')
        util.pickle_dump(jump_t,'jump_t.pkl')
    pyplot.plot(jump_t,jump_d,'','s',xlabel='time',ylabel='jump distance',save='msd_jump')
    hist_s,xs,max_hist=histogram(jump_d,bins=20)
    pyplot.plot(xs,hist_s,xlabel='distance',ylabel='count',save='jump_hist')
def msd_jump(VW,L):
    reload(diffuse)
    fid = open('large.xyz','r')
    M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
    crystal = M.cord_auto(['V','W'])
    bcc = np.zeros((777,432,1))
    x = range(50,500,50)
    delta = 50
    for frame in range(bcc.shape[0]):
        for i in range(bcc.shape[1]):
            if crystal[frame][i][0] > L[0]:
                bcc[frame][i]=0
            else:
                bcc[frame][i]=1
    jumps = []
    for i in x:
        jumps.extend(diffuse.crystal_jump(VW,bcc,L,i,delta))
    print len(jumps)
    import MD.plot.pyplot_eps as pyplot_eps
    util.pickle_dump(jumps,'jump.pkl')
    hist_s,xs=histogram_reg(jumps,bins=20)
    pyplot_eps.plot_bar(xs,hist_s,xlabel=r'$\sigma$',ylabel='count',save='jump_hist')
############################################
def msd_index(VW,L):
    k=100
    time = 500

    surface =  util.pickle_load('surface.pkl')
    DV =  util.pickle_load('DV_s.pkl')
    DW =  util.pickle_load('DW_s.pkl')
    CV =  util.pickle_load('CV.pkl')
    CW =  util.pickle_load('CW.pkl')

    defects = []
    print DV.shape
    print DW.shape
    for i in range(DV[k].shape[0]):
        if DV[k][i] == -1:
            index, d = close_neighbors_point(VW[k],CV[i],L)
            defects.append(index)
    for i in range(DW[k].shape[0]):
        if DW[k][i] == -1:
            index, d = close_neighbors_point(VW[k],CW[i],L)
            defects.append(index)
    print len(defects)


    fid = open('large.xyz','r')
    M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
    crystal = M.cord_auto(['V','W'])
    bcc = np.zeros((777,432,1))
    for frame in range(bcc.shape[0]):
        for i in range(bcc.shape[1]):
            if crystal[frame][i][0] > L[0]:
                bcc[frame][i]=0
            else:
                bcc[frame][i]=1
    s = []
    for i in range(len(bcc[k])):
            if surface[k][i]:
                s.append(i)
    gel = []
    for i in range(len(bcc[k])):
            if bcc[k][i] ==0 and surface[k][i] == 0:
                gel.append(i)
    crystal = []
    for i in range(len(bcc[k])):
            if bcc[k][i] == 1:
                crystal.append(i)

    #diffuse.index_msd(VW,s,L,time,save='msd_index_%i_surface_%i'%(k,len(s)))
    x,msd=diffuse.index_msd(VW,gel,L,time,save='msd_index_%i_gel_%i'%(k,len(gel)))
    x2,msd2=diffuse.index_msd(VW,crystal,L,time,save='msd_index_%i_solid_%i'%(k,len(crystal)))
    util.pickle_dump([x,msd,x2,msd2],'msd_index.pkl')
    #diffuse.index_msd(VW,defects,L,time,save='msd_index_%i_defects_%i'%(k,len(defects)))
    #r = list(set(s)-set(defects))
    #diffuse.index_msd(VW,r,L,time,save='msd_index_%i_nodefect%i'%(k,len(r)))
    #print len(defects)
    #print len(s)
    import MD.plot.pyplot_eps as pyplot_eps
    pyplot_eps.plot2(x,msd,x2,msd2,xlabel='t',ylabel='msd',label1='gel',label2='solid',save='msdgelsolid')
############################################
def jump_lattice(VW,L):
    reload(diffuse)
    k=100
    time = 500

    surface =  util.pickle_load('surface.pkl')
    CV =  util.pickle_load('CV.pkl')
    CW =  util.pickle_load('CW.pkl')
    DV =  util.pickle_load('DV_s.pkl')
    DW =  util.pickle_load('DW_s.pkl')
    defects = []
    for i in range(DV[k].shape[0]):
        if DV[k][i] == -1:
            index, d = close_neighbors_point(VW[k],CV[i],L)
            defects.append(index)
    for i in range(DW[k].shape[0]):
        if DW[k][i] == -1:
            index, d = close_neighbors_point(VW[k],CW[i],L)
            defects.append(index)

    proj = np.zeros((CV.shape[0]*2,CV.shape[1]))
    for i in range(CV.shape[0]):
        proj[i] = CV[i]
    for i in range(CW.shape[0]):
        proj[i+CV.shape[0]] = CW[i]


    fid = open('large.xyz','r')
    M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
    crystal = M.cord_auto(['V','W'])
    bcc = np.zeros((777,432,1))
    for frame in range(bcc.shape[0]):
        for i in range(bcc.shape[1]):
            if crystal[frame][i][0] > L[0]:
                bcc[frame][i]=0
            else:
                bcc[frame][i]=1
    s = []
    for i in range(len(bcc[k])):
            if surface[k][i]:
                s.append(i)
    gel = []
    for i in range(len(bcc[k])):
            if bcc[k][i] ==0 and surface[k][i] == 0:
                gel.append(i)
    crystal = []
    for i in range(len(bcc[k])):
            if bcc[k][i] == 1:
                crystal.append(i)

    diffuse.projection_jump(VW[k],s,proj,L,time,save='msd_index_%i_surface_%i'%(k,len(s)))
    diffuse.projection_jump(VW[k][:216],gel[:216],CW,L,time,save='msd_index_%i_gelV_%i'%(k,len(gel)))
    diffuse.projection_jump(VW[k][216:],gel[216:],CV,L,time,save='msd_index_%i_gelW_%i'%(k,len(gel)))
    diffuse.projection_jump(VW[k],crystal,proj,L,time,save='msd_index_%i_solid_%i'%(k,len(crystal)))
    diffuse.projection_jump(VW[k][:216],crystal[:216],CW,L,time,save='msd_index_%i_solidW_%i'%(k,len(crystal)))
    for i in range(len(crystal)):
        crystal[i] = crystal[i] - CV.shape[0]
    diffuse.projection_jump(VW[k][216:],crystal[216:],CV,L,time,save='msd_index_%i_solidV_%i'%(k,len(crystal)))
    diffuse.projection_jump(VW[k],defects,proj,L,time,save='msd_index_%i_defects_%i'%(k,len(defects)))
    r = list(set(s)-set(defects))
    diffuse.projection_jump(VW,r,proj,L,time,save='msd_index_%i_nodefect%i'%(k,len(r)))
    print len(defects)
    print len(s)
############################################
## 	g(r)
############################################
def distance_distribution(V,W,L,n_frames=15,n_start=1):
    #finds end to end distances and plots a histogram with that data
    def hist(M1,M2,save_name,label):
        distance=particle_distance(M1,M2,L)
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
        try:
            print 'Frame',i,'Temp',temp[i]
        except:
            print i
        x=np.arange(i,step_range+i,delta)
        #Find he number of connections at specific points x
        remain = life.lifetime(C[x],G[x],L)
        x=np.arange(i,step_range+i,delta)*50000
        print remain
        pyplot_eps.plot(x,remain,xlabel='time', ylabel='remaining connections',
                save='lifetime%i'%i)
        plt.close()
######################
# FInd the solid particles in a simulation
#####################
def solid_particles_vmd(VW,V,W,L,M,skip=15,step=5e4,rcut=False):
    VW_names = M.get_names(['V','W'])
    reload(sp)
    sp.solid_particles_VW(VW,V,W,L,VW_names,skip=skip,step=5e4,rcut=False)
    #sp.solid_particles_surround(VW,V,W,L,VW_names,skip=skip,step=5e4,rcut=False)
    #if os.path.exists('bcccrystals.xyz') == False:
    #    sp.solid_particles(VW,V,W,L,skip=skip,step=5e4,rcut=False)
    #sp.crystal_size(L,VW.shape[0]/skip,VW.shape[1])
#####################################################
#find number of a specified pair of particles that are within r_cut of each
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
#####################################################
# Find Static Structure Facto
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor(VW, L, n_frames=10, n_start=1, 
        filters=0.05, dh=0.05,save='sf'):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    stmp,qx = sf.sfactor(VW[n_start:n_start+n_frames],L=L,l=20) 
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
    #Plot the graph
    xlabel = '$|\\vec{q}$ $|$'
    ylabel = '$S(\\vec{q}$ $)$ '
    pyplot.plot(Q, S, linestyle='', marker='x', 
            xlabel=xlabel, ylabel=ylabel, save=save)
    util.pickle_dump(Q,'Vsfx%i.pkl'%n_start)
    util.pickle_dump(S,'Vsfy%i.pkl'%n_start)
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
    Q, S, primitive_vectors = sf.sf_filter(Q, S, primitive_vectors, filters=filters)
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
    rcut = 22
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
            label2='other', save=('mixing_avg_cut%.1f'%rcut), showleg=True)
#####################################################
# Look at mixing of species of spheres
#####################################################
def mixing_frame(V,W,L,frame,rcut=False):
    #The total number of frames we are going to look at
    #rcut = 19
    rcut = 22
    k = frame
    Same = []
    Other = []
    print k
    for k in range(frame,frame+10):
        for i in range(V[k].shape[0]):
            Same.append(len(nearest_neighbors_point(V[k],V[k][i],L,rcut)[0]))
            Same.append(len(nearest_neighbors_point(W[k],W[k][i],L,rcut)[0]))
            Other.append(len(nearest_neighbors_point(W[k],V[k][i],L,rcut)[0]))
            Other.append(len(nearest_neighbors_point(V[k],W[k][i],L,rcut)[0]))
    hist_same, x_s = histogram_reg(Same,8)
    hist_other, x_o = histogram_reg(Other,8)
    pyplot.plot_bar(x_o, hist_other, label='other', save=('mix_other%i'%k), showleg=True)
    pyplot.plot_bar(x_s, hist_same, label='Same',
            save=(('mix_same_frame%i_rcut%.1f')%(k,rcut)), showleg=True,color='r')
## \brief read mylog file and makes a plot of the energy
#
# \returns x coordinates of time step
# \returns y whatever value of the row you are looking for
#
# \param row the row you would like to grab data from 
def mylog(row = 2):
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
def run_defect():
    #setup
    #print out directory
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    try:
        var = MD.util.get_path_variables()
    except:
        var = {'ndna':60}
    print var
    M=MD.ReadCord()
    L = M.box_length
    last = M.frames
    try:
        V = util.pickle_load('V.pkl')
        W = util.pickle_load('W.pkl')
        R = util.pickle_load('R.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
        V=M.cord_auto(['V'])
        W=M.cord_auto(['W'])
        R=M.cord_auto(['R'])
        VW=M.cord_auto(['V','W'])
        VW,V,W = drift_remove_all(VW,V,W,L)
        R = drift_remove(R,L)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(W,'W.pkl')
        util.pickle_dump(VW,'VW.pkl')
        util.pickle_dump(R,'R.pkl')
    x = range(0,last,last/5)
    #mylog()
    #plt.close()
    solid_particles_vmd(VW,V,W,L,var,skip=100)
    plt.close()
    #species_mixing(V,W,L,n_finish=last)
    #plt.close()
    #species_mixing_fcc(VW,L,n_finish=last)
    #plt.close()
    for i in x:
        #print "finding s(q)"
        structure_factor(VW, L, var, n_start=i,save='sf'+str(i))
        print "finding g(s)"
        distance_distribution(V,W,L,n_start=i)
        plt.close()
    msd(VW,L)
    plt.close()
    msd(R,L,save='msd_R')
    plt.close()
    #end_end(V,W,M,L,var)
    #find_networks(M, VW, L, var, last,delta=25)
    #plt.close()
    #sphere_rotations(M,W,L,'N')
def run_test():
    #print out directory dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    print L
    last = M.frames
    delta = 25
    try:
        V = util.pickle_load('V.pkl')
        W = util.pickle_load('W.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
        V=M.cord_range(['V'],delta=delta)
        W=M.cord_range(['W'],delta=delta)
        VW=M.cord_range(['V','W'],delta=delta)
        #V=M.cord_auto(['V'])
        #W=M.cord_auto(['W'])
        #VW=M.cord_auto(['V','W'])
        VW_index=M.get_names(['V','W'])
        VW,V,W = drift_remove_all(VW,V,W,L,VW_index)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(W,'W.pkl')
        util.pickle_dump(VW,'VW.pkl')
    if V.shape[0]>5:
        x = range(0,V.shape[0],V.shape[0]/3)
    else:
        x = range(V.shape[0])
    delta = 5
    #step,temp = mylog(row = 2)
    #plt.close()
    msd(VW,L,time_scale=delta)
    print V
    #msd_phases(VW,L)
    #msd_jump(VW,L)
    #sd_index(VW,L)
    #jump_lattice(VW,L)
    #solid_particles_vmd(VW,V,W,L,M,skip=25)
    #solid_particles_vmd(VW[1000:],V[1000:],W[1000:],L,var,skip=10)
    #plt.close()
    #species_mixing(V,W,L,n_finish=last)
    #plt.close()
    #end_end(V,W,M,L)
    ##plt.close()
    for i in x:
        print "finding s(q)"
        structure_factor(VW, L, n_start=i,save='sf'+str(i))
        #structure_factor(V, L, n_start=i,save='sfV'+str(i))
        #structure_factor(W, L, n_start=i,save='sfW'+str(i))
        print "finding g(s)"
        distance_distribution(V,W,L,n_start=i)
        plt.close()
        #mixing_frame(V,W,L,i)
        #plt.close()
    #msd(VW,L)
    #plt.close()
    #temp=1.195
    #find_lifetime(M,L,x,temp,step_range=30,delta=5)
    #plt.close()
    #find_networks(M, VW, L, 100,delta = 49)
    #plt.close()
    #sphere_rotations(M,W,L,'N')
def run_blah():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    print L
    last = M.frames
    f=150
    def distance_distribution(V,W,L,n_start=1,save=''):
        #finds end to end distances and plots a histogram with that data
        def hist(M1,M2,save_name,label):
            distance=particle_distance(M1,M2,L)
            dist = []
            for d in distance:
                if d>1 and d<25:
                    dist.append(d)
            hist_s,xs=histogram_reg(dist,bins=25)
            return xs,hist_s
        ########################################
        #Plot Data for AA,AB,BB distances
        save_name=save+'_time_'+'%i'%(n_start)
        max_hist=0
        AB=hist(V,W,save_name,'A-B')
        AA=hist(W,W,save_name,'A-A')
        BB=hist(V,V,save_name,'B-B')
        pyplot_eps.plot3(AA[0],AA[1],AB[0],AB[1],BB[0],BB[1],'s','g(s)',
                'A-A','A-B','B-B',save='nnplot_'+save_name,showleg=True)
    V = util.pickle_load('V.pkl')
    W = util.pickle_load('W.pkl')
    VW = util.pickle_load('VW.pkl')
    surface =  util.pickle_load('surface.pkl')
    fid = open('large.xyz','r')
    M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
    crystal = M.cord_auto(['V','W'])
    bcc = np.zeros((777,432,1))
    for frame in range(bcc.shape[0]):
        for i in range(bcc.shape[1]):
            if crystal[frame][i][0] > L[0]:
                bcc[frame][i]=0
            else:
                bcc[frame][i]=1
    s = []
    for i in range(len(surface[f])):
            if surface[f][i]:
                s.append(i)
    ###############################
    gel = []
    for i in range(len(bcc[f])):
            if bcc[f][i] ==0 and surface[f][i] == 0:
                gel.append(i)
    g = np.zeros((15,len(gel),3))
    gV = np.zeros((15,len(gel),3))
    gW = np.zeros((15,len(gel),3))
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            if gel[i] < V.shape[1]:
                gV[k][i] = V[k+f][gel[i]]
            else:
                gW[k][i] = W[k+f][gel[i]-V.shape[1]]
    distance_distribution(gV,gW,L,f,'gel')
    print 'gel finished'
    ###########################
    crystal = []
    for i in range(len(bcc[f])):
        if bcc[f][i] == 1:
            crystal.append(i)
    g = np.zeros((15,len(crystal),3))
    gV = np.zeros((15,len(crystal),3))
    gW = np.zeros((15,len(crystal),3))
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            if crystal[i] < V.shape[1]:
                gV[k][i] = V[k+f][crystal[i]]
            else:
                gW[k][i] = W[k+f][crystal[i]-V.shape[1]]
    distance_distribution(gV,gW,L,f,'crystal')
    print 'crystal finished'
def run_blah2():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    surface =  util.pickle_load('surface.pkl')
    fid = open('large.xyz','r')
    M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
    crystal = M.cord_auto(['V','W'])
    ###########################
    bcc = np.zeros((777,432,1))
    for k in range(bcc.shape[0]):
        for i in range(bcc.shape[1]):
            if crystal[k][i][0] > L[0]:
                bcc[k][i]=0
            else:
                bcc[k][i]=1
    crystal = []
    for k in range(len(bcc)):
        c=[]
        for i in range(len(bcc[k])):
            if bcc[k][i] == 1:
                c.append(i)
        crystal.append(len(c))
    s = []
    for k in range(len(surface)):
        surf = []
        for i in range(len(surface[k])):
            if surface[k][i]:
                surf.append(i)
        s.append(surf)
    util.pickle_dump(s,'surface_diffusion.pkl')
    util.pickle_dump(crystal,'crystal_diffusion.pkl')
    lamda=4
    D=0.05
    ns = []
    dt = 10
    dmdt = []
    C = lamda**2/(D*24)
    start = 100
    finish = 350 
    for k in range(start,finish,dt):
        n = 0
        for i in range(dt):
            n += len(s[i+k]) 
        ns.append(n/dt)

    count = 0
    for i in range(start,finish,dt):
        dm = crystal[i]-crystal[i-dt]
        dmdt.append(C*dm/ns[count])
        count += 1
    x = np.arange(start,finish,dt)
    pyplot.plot(x,dmdt,save='diffusion')

#For multiple directories
if __name__ == '__main__':
    for f in sorted(os.listdir("./")):
        if os.path.isdir(f):
            os.chdir(f)
            run_test()
            #if os.path.isdir('cont'):
                         # os.chdir('cont')
     #           run_test()
    #            os.chdir('../')
            #run_defect()
            #try:
            #    run_defect()
            #except:
            #    print "Someting Failed"
            #print '##\nfinished with:',f
            os.chdir('../')
#for single directories
if __name__ == '__main__':
#    #run_debug()
#    #run_all()
    run_test()
    #run_blah()
