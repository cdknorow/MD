# -*- coding: utf-8 -*-
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
import MD.base.points as points
import MD.base.clusters as clusters
import MD.plot.pyplot as pyplot
import MD.plot.pyplot_eps as pyplot_eps
from MD.analysis.particle_distance import particle_distance
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import second_nearest_neighbors_index
from MD.analysis.nearest_neighbor import count_neighbors_index
from MD.analysis.nearest_neighbor import min_distance
from MD.plot.histogram import histogram
from MD.plot.histogram import histogram_normal
from MD.plot.histogram import histogram_reg
from MD.analysis.rotation import single_rotation_box ############ 
import MD.unit.make_bcc as mb
import draw_cube as dw
# Run from inside a folder
#########
reload(MD)
reload(pyplot) ########################################3
def dump_xyz(V,save='simple.xyz'):
    out = open(save,'w')
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
## returns the average position over a few frames
############################################
def Average(VW,L,n_start,n_frames,write=True):
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
        fid = open('star.xyz','w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()
    return Average
## \brief read the bakos_index file an assign cubes to a matrix
# \returns array containing matrix of sorted bakos cubes
#
# \param VW: VW index
# \param frame: frame to search
def filter_bakos(VW,frame,select=False):
    save_name='_time_'+'%i'%(frame)
    max_hist=0
    fid = open('animate/bakos_index%i.txt'%frame,'r')
    M = fid.readlines()
    A =[]
    B = []
    C = []
    D = []
    print len(M)
    print VW.shape
    index  =[[] for i in range(4)]
    if select != False:
        for i,line in enumerate(M):
            if int(line) in select:
                A.append(VW[i])
            else:
                B.append(VW[i])
        A = np.array(A)
        B = np.array(B)
        return [A,B]

    else:
        for i,line in enumerate(M):
            if int(line) == 0:
                A.append(VW[i])
                index[0].append(i)
            if int(line) == 1:
                B.append(VW[i])
                index[1].append(i)
            if int(line) == 2:
                C.append(VW[i])
                index[2].append(i)
            if int(line) == 3:
                D.append(VW[i])
                index[3].append(i)
        A = np.array(A)
        B = np.array(B)
        C = np.array(C)
        D = np.array(D)
        return [A,B,C,D],index
## \brief read the bakos_index file retrieving color indexes
# \returns array containing color index of cube
#
# \param frame: frame to search

def read_color_index(frame = ''):
    if frame == '':
        fid = open('animate/bakos_index.txt','r')
    else:
        fid = open('animate/bakos_index%i.txt'%frame,'r')
    M = fid.readlines()
    A = []
    for line in M:
        A.append(float(line))
    return A
## \brief read mylog file and makes a plot of the energy
#
# \returns x coordinates of time step
# \returns y whatever value of the row you are looking for
#
# \param row the row you would like to grab data from 
def mylog_volume(L,row = 2):
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
    x = []
    for i in L:
        x.append(i[0])
    x[-1] = 35
    #make a simple plot
    label = re.sub("_"," ",text[0].split()[row])
    save = text[0].split()[row]
    print x
    pyplot.plot(x, y, xlabel='Volume', ylabel=label, save=save)
    fid = open('box_size','w')
    for i in range(len(x)):
        print x[i]
        print y[i]
        s = '%.2f %s \n'%(x[i],y[i])
        fid.write(s)
    return x, y
############################################
## 	g(r)
############################################
def distance_distribution_bakos(V,L,x,smooth=1):
    #finds end to end distances and plots a histogram with that data
    def bakos_dist(M,N,L,rcut = 30):
        distance =[]
        print L
        for i in range(M.shape[0]):
            for j in range(N.shape[0]):
                distance.append(points.dist(M[i],N[j],L)[0])
                if distance[-1] > rcut or distance[-1]<1:
                    del distance[-1]
        return distance
    def hist(distance):
        hist_s,xs,max_hist=histogram_normal(np.array(distance),bins=np.arange(5,30,0.05))
        #hist_s,xs,max_hist=histogram_normal(np.array(distance),bins=30)
        #normalize the function with respect to an ideal gas
        return xs,hist_s
    ########################################
    #Plot Data for AA,AB,BB distances
    AA = []
    AB = []
    AC = []
    AD = []
    for frame in x:
        save_name='_time_'+'%i'%(frame)
        max_hist=0
        fid = open('animate/bakos_index%i.txt'%frame,'r')
        M = fid.readlines()
        A =[]
        B = []
        C = []
        D = []
        for i,line in enumerate(M):
            if int(line) == 1:
                A.append(V[frame][i])
            if int(line) == 2:
                B.append(V[frame][i])
            if int(line) == 3:
                C.append(V[frame][i])
            if int(line) == 0:
                D.append(V[frame][i])

        A = np.array(A)
        B = np.array(B)
        C = np.array(C)
        D = np.array(D)

        AA.extend(bakos_dist(A,A,L[frame]))
        AB.extend(bakos_dist(A,B,L[frame]))
        AC.extend(bakos_dist(A,C,L[frame]))
        AD.extend(bakos_dist(A,D,L[frame]))
    print "neigbors for AA"
    print len(AA)
    print "neigbors for AB"
    print len(AB)
    print "neigbors for AC"
    print len(AC)
    print "neigbors for AD"
    print len(AD)
    AAtotal = hist(AA)
    ABtotal = hist(AB)
    ACtotal = hist(AC)
    ADtotal = hist(AD)
    fid = open('nnbakos.txt','w')
    for i in range(AAtotal[0].shape[0]):
        fid.write('%.3f %.6f %.6f %.6f %.6f\n'%
                (AAtotal[0][i],AAtotal[1][i],ABtotal[1][i],ACtotal[1][i],ADtotal[1][i]))
    fid.close()



    pyplot.plot4(AAtotal[0],AAtotal[1],ABtotal[0],ABtotal[1],ACtotal[0],ACtotal[1],ADtotal[0],ADtotal[1],xlabel='s',ylabel='g(s)',save='nnplot_'+save_name,
            label1='A-A',label2='A-B',label3='A-C',label4='A-D',showleg=True)
    #pyplot.plot(AB[0],AB[1],xlabel='s',ylabel='g(s)',save='nnplot_'+save_name)
    #pyplot.plot(AC[0],AC[1],xlabel='s',ylabel='g(s)',save='nnplot_'+save_name)
    #pyplot.plot(AD[0],AD[1],xlabel='s',ylabel='g(s)',save='nnplot_'+save_name)
#####################################################
# Find Static Structure Facto
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor_bakos(VW, L, n_frames=10, n_start=1, 
        filters=0.05, dh=0.05,save='sf',l=10):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    A = Average(VW,L,n_start,n_frames)
    Bakos = filter_bakos(A,n_start)
    for count,A in enumerate(Bakos):
        dump_xyz(np.array([A]),('simple%i.xyz'%count))
        stmp,qx = sf.sfactor(np.array([A]),L=L[n_start],l=l)
        S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
        #Plot the graph
        xlabel = '$|\\vec{q}$ $|$'
        ylabel = '$S(\\vec{q}$ $)$ '
        #sp4 = 0.634698
        #sp12 = 0.530493
        #sc 1.121997
        pyplot_eps.plot_sf(Q, S,0.817082, linestyle='', marker='x', 
                xlabel=xlabel, ylabel=ylabel, save=save+'_%i'%count,xlen=20)
        util.pickle_dump(Q,'Vsfx%i.pkl'%n_start)
        util.pickle_dump(S,'Vsfy%i.pkl'%n_start)
        print 'finished reconstruction'
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
                filters=filters,save=save+'_%i'%count)
        import MD.analysis.reconstruct_lattice as rl
        import MD.unit.make_bcc as mb
        con = True
        while con:
            a1,a2,a3 =  rl.reconstruct_lattice(Q,primitive_vectors,save+'recipricol%i.txt'%count)
            print a1,a2,a3
            mb.make_bcc(a1,a2,a3,A[0],L[n_start])
            if raw_input('Reselect? (y)/N') == 'N':
                con = False
    return Q, S, primitive_vectors
#\brief find the gauss map of neighbors surrounding NC
# but do not allighn with axis of cube
def neighbor_no_gauss(VW, Z,L, frames, step=5e4,rcut=25,neighbors=8):
    def draw_cubes(A,A_Z,color_index,i,N,select=[2]):
        color = ['blue','red','green','orange','purple']
        if i in select:
            c = A[N]
            #mapping from A to VW
            v1 = A_Z[N*6] - c
            v2 = A_Z[N*6+2] - c
            v3 = A_Z[N*6+4] - c
            if points.dist_np(v1+c,A[N])[0]<7:
                if points.dist_np(v2+c,A[N])[0]<7:
                    if points.dist_np(v3+c,A[N])[0]<7:
                        c = d[1]
                        dw.draw_cube(fid,c,v1,v2,v3,color[int(color_index[N])])
            c = A[i]
            #mapping from A to VW
            v1 = A_Z[i*6] - c
            v2 = A_Z[i*6+2] - c
            v3 = A_Z[i*6+4] - c
            if points.dist_np(v1+c,A[i])[0]<7:
                if points.dist_np(v2+c,A[i])[0]<7:
                    if points.dist_np(v3+c,A[i])[0]<7:
                        dw.draw_cube(fid,np.array([0,0,0]),v1,v2,v3,color[int(color_index[i])])
    from MD.analysis.nearest_neighbor import count_neighbors_index
    for k in frames:
        #Find Structure Factor
        print 'finding average'
        A = Average(VW,L,k,5)
        A_Z = Average(Z,L,k,5,write=False)
        print 'finished average'
        #####################################
        gauss_map = []
        N_location = []
        fid = open('animate/gauss_cube.tcl','w')
        color_index = read_color_index(frame=k)
        for i in range(len(A)):
            select = [11]
        #for i in select:
            #we must rotate about a specific cube reference frame
            for N in count_neighbors_index(A,i,L[k],count=neighbors,rcut=rcut)[0]:
                d = points.dist(A[i],A[N],L[k])
                N_location.append(d[1])
                draw_cubes(A,A_Z,color_index,i,N,select=select)
        fid.close()
        gauss_map.append(N_location)
        #########
        fid = open(('gaussmap_frame_%i_neighbors_%i_all.xyz'%(k,neighbors)),'w')
        max_gauss = 0
        fid.write('%i\n%i\n'%(len(gauss_map[0])+len(gauss_map)+4-1,k))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        for i in gauss_map[0]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        fid.close()
        #####################################
        Bakos,Bindex = filter_bakos(A,k)
        for count,A in enumerate(Bakos):
            print A[1]
            dump_xyz(np.array([A]),('simple%i.xyz'%count))
            gauss_map = []
            N_location = []
            for i in range(len(A)):
                #we must rotate about a specific cube reference frame
                for N in count_neighbors_index(A,i,L[k],count=neighbors,rcut=rcut)[0]:
                    d = points.dist(A[i],A[N],L[k])
                    N_location.append(d[1])
            gauss_map.append(N_location)
            #find the clusters of particles 
            reload(clusters)
            CM = clusters.cluster(gauss_map[0],rcut=2,cluster_cut=4)
            #########
            fid = open(('gaussmap_frame_%i_neighbors_%i_index%i.xyz'%(k,neighbors,count)),'w')
            max_gauss = 0
            fid.write('%i\n%i\n'%(len(gauss_map[0])+4+len(CM),k))
            fid.write('E  1 0 0\n')
            fid.write('E  0 1 0\n')
            fid.write('E  0 0 1\n')
            fid.write('V  0 0 0\n')
            #write out center of mass of clusters
            for i in CM:
                fid.write('R  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            #write out location of all particles
            for i in gauss_map[0]:
                fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            fid.close()
            #recipritocl lattice
            #a1 = CM[5]
            #a2 = CM[9]
            #a3 = CM[4]
            #sp4/sim7/cont/cont
            #a1 = np.array([10.722,2.017,-0.805])
            #a2 = np.array([-2.32,2.619,-11.189])
            #a3 = np.array([2.652,10.031,-4.031])
            #/sp4/sim3/cont
            a1 = np.array([-1.732,12.1,-1.685])
            a2 = np.array([-8.531,3.347,-8.682])
            a3 = np.array([-10.825,4.342,3.47])
            print A[1]
            mb.make_bcc(a1,a2,a3,A[1],L[k],name='recipricol_index%i.xyz'%count)
#\brief find the gauss map of neighbors surrounding NC
# but do not allighn with axis of cube
def neighbor_correlations(VW, Z,L, frames, step=5e4,rcut=25,neighbors=8):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    correlations = []
    for k in frames:
        #Find Structure Factor
        A = Average(VW,L,k,5)
        #####################################
        color = read_color_index(frame=k)
        gauss_map = []
        N_location = []
        count_color=0
        count_total=0
        for i in range(len(A)):
            #we must rotate about a specific cube reference frame
            for N in count_neighbors_index(A,i,L[k],count=neighbors,rcut=rcut)[0]:
                if color[N] == color[i]:
                    count_color+=1
                count_total+=1
        correlations.append(count_color/float(count_total))
    pyplot.plot(frames,correlations,xlabel='frames',ylabel='c(nn)',save='correlations')

## \brief find the end to end distance of polymers in f-star systems
# \returns average stretching of polymer
#
# \param M: cord reader
def end_end(M,L):
    m=M.cord_auto(['M'])
    edge = [ 0.0  for i in range(m.shape[0])]
    for k in range(0,m.shape[0]):
        for i in range(0,m.shape[1],2):
            edge[k] += points.dist(m[k][i],m[k][i+1],L[k])[0]
        edge[k]/=(m.shape[1]/2)
    pyplot.plot(range(m.shape[0]),edge,save='polymer_length')
    return edge


##\brief :wq
#################################################
# What to run
#################################################
def run_():
    #setup
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] =L_cont[-2]
    last = M.frames
    try:
        V = util.pickle_load('V.pkl')
        Z = util.pickle_load('Z.pkl')
    except:
        V=M.cord_auto(['V','W'])
        Z=M.cord_auto(['Z'])
        VW_index=M.get_names(['V','W'])
        Z_index=M.get_names(['Z'])
        #V = drift_remove(V,L,VW_index)
        #Z = drift_remove(Z,L,Z_index)
        #dump_xyz(V)
        #dump_xyz(Z)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(Z,'Z.pkl')
    #end_end(M,L_cont)
    #cubics(M,VW,L,var['ndna'])
    #rotations(V,Z,L_cont,'Z')
    #plt.close()
    #msd(V,L)
    #plt.close()
    #x = range(0,V.shape[0],delta)
    #solid_particles_bcc(V,L_cont,x)
    #solid_particles_sc(V,L_cont,x)
    ##cube_angle(Z,V,L_cont)
    #x = range(0,V.shape[0],delta)
    #x = range(4500,4550)
    #distance_distribution_bakos(V,L_cont,x,smooth=1)
    #x = [4500]
    x = [130]
    neighbor_no_gauss(V, Z,L_cont, x, rcut=20,neighbors=12)
        #for i in x:
    #    structure_factor_bakos(V, L_cont, n_start=i,save='sf'+str(i),
    #            n_frames = 5,l=25, filters=0.05)
    #    plt.close()
        #print "finding g(s)"
        #distance_distribution(V,L_cont[i],n_start=i)
        #plt.close()
    #x = range(0,90,10)
    #neighbor_correlations(V, Z,L_cont, x,rcut=22,neighbors=12)
#For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            run_()
#            print '##\nfinished with:',f
#            os.chdir('../')
#For single directories
if __name__ == '__main__':
    run_()

