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
from MD.analysis.nearest_neighbor import count_neighbors_distance
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
########################################3
def dump_xyz(V):
    out = open('simple.xyz','w')
    print V.shape
    for i in range(V.shape[0]):
        out.write(('%i\n\n')%(V.shape[1]))
        for j in range(V.shape[1]):
            out.write(('V %.2f %.2f %.2f\n')%(V[i][j][0],V[i][j][1],V[i][j][2]))
############################################
## returns the average position over a few frames
############################################
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
## msd
############################################
def msd_compress(VW,L,frame,step=1):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x_total = []
    msd_total = []
    for i in range(1,len(frame)):
        print frame[i],L[frame[i-1]]
        if frame[i]-frame[i-1] >2:
            x,m=msd(VW[frame[i-1]:frame[i]],L[frame[i-1]],step=step)
            x_total.extend(x+frame[i-1])
            msd_total.extend(m**0.5)
    print x_total
    pyplot.plot(x_total,msd_total,xlabel='time',ylabel='msd',save='MSDtime')
    util.pickle_dump(msd_total,'msd.pkl')
def print_box_volume(L_cont,delta=10):
    fid = open('box_length.txt','w')
    fid.write('frame  L\n') 
    for i in range(0,len(L_cont)):
        fid.write(('%i %.2f %.2f %.2f\n')%(i,L_cont[i][0],L_cont[i][1],L_cont[i][2]))
    fid.close()
#\brief find the polymer stretching and center to end distances
def polymer_inter_gauss(M,VW, L, frames, rcut=12, step=5e4):
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
    min_D_same = rcut
    print P.shape
    Lstart = 0
    for k in frames:
        Lstart = L[k][0]
        location = []
        location_other = []
        print k
        try:
            for i in range(VW.shape[1]):
                #we must rotate about a specific cube reference frame
                if i in range(1):
                    V = VW[k][i]
                    same_set = []
                    for j in range(1,ndna,2):
                        location.append(points.dist(V,P[k][j+i*ndna],L[k])[1])
                        d = points.dist(V,P[k][j+i*ndna],L[k])[0]
                        if d < min_D_same:
                            min_D_same = d
                        same_set.append(j+i*ndna)
                    for poly in range(1,P.shape[1],2):
                        print poly
                        if poly not in same_set:
                            d = points.dist(V,P[k][poly],L[k])[0]
                            if d < rcut:
                                if d< min_D:
                                    min_D = d
                                location_other.append(points.dist(V,P[k][poly],L[k])[1])
        except:
            print 'ERROR: Failed to rotate frame!'
        gauss_map.append(location)
        gauss_inter.append(location_other)
    #fid.close()
    ###########
    fid = open('gaussmap_polymers_total%.2f.xyz'%Lstart,'w')
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
    delta = 6/4.0 
    for k in range(len(gauss_map)):
        for i in gauss_map[k]:
            d = points.dist(np.array([0,0,0]),i,L[frames[k]])[0]
            if d < min_D_same+delta:
                fid.write('A  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            elif d < min_D_same+delta*2:
                fid.write('B  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            elif d < min_D_same+delta*3:
                fid.write('C  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            else:
                fid.write('D  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
            count+=1
    fid.close()
    fid = open('gaussmap_polymers_inter%.2f.xyz'%Lstart,'w')
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
#\brief find the angle the polymer makes with the grafting surface
#
# uses law of cosines with
#           a = polymer end to end
#           b = center to start of polymer
#           c = polymer center to end
def polymer_angle(M,VW, L, frames, rcut=12, step=5e4):
    try:
        P = util.pickle_load('M.pkl')
    except:
         P = M.cord_auto(['M'])
         util.pickle_dump(P,'M.pkl')
    ndna = P.shape[1]/VW.shape[1]
    frames_theta = []
    for k in frames:
        print k
        theta = []
        try:
            for i in range(VW.shape[1]):
                #we must rotate about a specific cube reference frame
                if i in range(VW.shape[1]):
                    V = VW[k][i]
                    for j in range(1,ndna,2):
                        a = points.dist(V,P[k][j+i*ndna-1],L[k])[0]
                        b = points.dist(P[k][j+i*ndna],P[k][j+i*ndna-1],L[k])[0]
                        c = points.dist(V,P[k][j+i*ndna],L[k])[0]
                        theta.append(np.arccos((a**2 + b**2 - c**2) / (2 * a * b))-math.pi/2.)
        except:
            print 'ERROR: Failed to rotate frame!'
        frames_theta.append(theta)
    multi_x = []
    multi_y = []
    label = []
    for i in range(len(frames)):
        hist, x, max_z = histogram_normal(frames_theta[i],np.arange(-1.0,2.0,0.2))
        print x
        print hist
        #pyplot.plot_bar(x,hist,xlabel='theta',xmax=2.0,ymax=0.3,
        #        ylabel='count',save='theta%.1f'%L[frames[i]][0])
        #plt.close()
        new_x = [0+x[0]/2.]
        for j in range(1,len(x)):
                new_x.append((x[j]+x[j])/2.0)
        print len(new_x)
        print len(hist)
        multi_x.append(new_x)
        multi_y.append(hist)
        label.append('L=%.1f'%L[frames[i]][0])
    fid = open('theta_hist.dat','w')
    fid.write('#theta ')
    for i in range(len(frames)):
        label.append('L=%.1f'%L[frames[i]][0])
        fid.write('  '+label[-1])
    fid.write('\n')
    for i in range(len(multi_x[0])):
        fid.write("  %.2f"%multi_x[0][i]) 
        for j in range(len(multi_y)):
            fid.write("  %.2f"%multi_y[j][i]) 
        fid.write('\n')
    pyplot.plot_multi(multi_x,multi_y,label,xlabel='theta',ylabel='count',
            title='Polymer theta', save='theta_hist')
#\brief draw the polymer for certain frames and NP
#
# NP = nanoparticles to use
# frames = frames to use
def draw_polymer(M,L, frames,NP=[1,2], rcut=1.0, step=5e4):
    print "getting cors" 
    P = M.cord_frames(['M'],frames)
    VW = M.cord_frames(['V','W'],frames)
    print "finished getting cords"
    ndna = P.shape[1]/VW.shape[1]
    #distance from center to edge of cube
    #draw edges of the cube
    def write_line(fid,p1,p2):
        s1 = "{%.2f %.2f %.2f}"%(p1[0],p1[1],p1[2])
        s2 = "{%.2f %.2f %.2f}"%(p2[0],p2[1],p2[2])
        fid.write(("draw arrow %s %s\n")%(s1,s2))
    c = ['green','blue','red']
    for count,k in enumerate(frames):
        print k
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
#\brief find the average center to end and end to end distance of polymers
#           relate these to the lattice size
def polymer_end_vs_lattice(M,VW, L, frames, rcut=1.0, step=5e4):
    try:
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
    e2e = []
    c2e = []
    lattice = []
    for k in frames:
        L_frame.append(k)
        location = []
        print k
        D = 0
        d_c2e = 0
        d_e2e = 0
        count = 0
        count_D = 0
        for i in range(VW.shape[1]):
            if i in range(10):
                D += sum(count_neighbors_distance(VW[k],i,L[k],rcut=30))/8.0
                count_D += 1
                for j in range(1,ndna,2):
                    d_c2e += points.dist(VW[k][i],P[k][j+i*ndna],L[k])[0]+0.5
                    d_e2e += points.dist(P[k][j+i*ndna-1],P[k][j+i*ndna],L[k])[0]+1.0
                    count += 1
        c2e.append(d_c2e/count)
        e2e.append(d_e2e/count)
        lattice.append((D/count_D))
    fid = open('e2e_lattice.dat','w')
    fid.write('#lattice e2e\n')
    for i in range(len(e2e)):
            fid.write('%.2f %.2f\n'%(lattice[i],e2e[i]))
    pyplot.plot(lattice,e2e,xlabel='a/2',
            ylabel='polymer stretching',save='e2e_lattice')
            #limx=(6,11),limy=(4,7))
    plt.close()
    pyplot.plot(lattice,c2e,xlabel='a/2-R',
            ylabel='center to poly',save='c2e_lattice')
    plt.close()
#################################################
#\brief find the neighbors of each polymer and plot them
# on a gauss map
def neighbor_no_gauss(VW, L, frames, step=5e4,rcut=25,count=8):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    gauss_map = []
    print VW.shape
    for k in frames:
        print k
        N_location = []
        for i in range(VW.shape[1]):
            counter = 0
            for N in count_neighbors_index(VW[k],i,L[k],count=count,rcut=rcut)[0]:
                counter+=1
                d = points.dist(VW[k][i],VW[k][N],L[k])[1]
                N_location.append(d)
            if counter>0:
                print counter
        gauss_map.append(N_location)
    #########
    fid = open(('gaussmap_no_neighbors%i.xyz'%count),'w')
    max_gauss = 0
    CM = clusters.cluster(gauss_map[0],rcut=2,cluster_cut=10)
    print CM
    for k in gauss_map:
        if len(k) > max_gauss:
            max_gauss = len(k)
    for k in range(len(gauss_map)):
        fid.write('%i\n%i\n'%(max_gauss+4+len(CM),frames[k]))
        fid.write('E  1 0 0\n')
        fid.write('E  0 1 0\n')
        fid.write('E  0 0 1\n')
        fid.write('V  0 0 0\n')
        #write out center of mass of clusters
        for i in CM:
            fid.write('R  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        for i in gauss_map[k]:
            fid.write('N  %.3f %.3f %.3f\n'%(i[0],i[1],i[2]))
        #for i in range(max_gauss - len(gauss_map[k]))):
        #    fid.write('N  0 0 0\n')
    fid.close()
#\brief find the neighbors of each polymer and plot them
# on a gauss map
def coordination_shell(VW, L, frames, step=5e4,rcut=25,sigma=1):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    shell = []
    print VW.shape
    for k in frames:
        print k
        N_location = []
        for i in range(VW.shape[1]):
            counter = 0
            #grab the distance to the second closest neighbor
            closest_vec = count_neighbors_index(VW[k],i,L[k],count=2,rcut=rcut)[1][1]
            np_cut = np.dot(closest_vec,closest_vec)**0.5
            N_location.append(len(nearest_neighbors_index(VW[k],i,L[k],rcut=np_cut+sigma)[0]))
            if counter>0:
                print counter
        shell.append(float(sum(N_location))/len(N_location))
    #########
    fid = open('cord_shell%.2f.dat'%sigma,'w')
    fid.write('#Frame shell L\n')
    for k, count in enumerate(shell):
        fid.write('%i %.2f %.2f\n'%(frames[k],count,L[frames[k]][0]))
    fid.close()
#\brief find the neighbors of each polymer and plot them
# on a gauss map
def coordination_gyration(VW, L, frames, step=5e4,rcut=25,count=8):
    from MD.analysis.nearest_neighbor import count_neighbors_index
    shell = []
    print VW.shape
    for k in frames:
        print k
        gyration = []
        for i in range(VW.shape[1]):
            counter = 0
            d = []
            for N in count_neighbors_index(VW[k],i,L[k],count=count,rcut=rcut)[0]:
                counter+=1
                d.append(points.dist(VW[k][i],VW[k][N],L[k])[0])
            Rg = 0
            average = sum(d)/len(d)
            for N in d:
                Rg += (N-average)**2
            gyration.append(Rg/len(d))
        shell.append(float(sum(gyration))/len(gyration)**0.5)
    #########
    fid = open('cord_gyrate%i.dat'%count,'w')
    fid.write('#Frame shell L\n')
    for k, count in enumerate(shell):
        fid.write('%i %.2f %.2f\n'%(frames[k],count,L[frames[k]][0]))
    fid.close()
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
        VW = util.pickle_load('VW.pkl')
    except:
         VW=M.cord_auto(['V','W'])
         util.pickle_dump(VW,'VW.pkl')
    #Find the msd of the system
    #VW = drift_remove(VW,L)
    dump_xyz(VW)
    x = []
    L_last = L_cont[0][0]
    count = 0
    for i,j in enumerate(L_cont):
        if j[0] != L_last:
            L_last = j[0]
            if i - 5 > 0:
                count +=1
                if count%3 == 0:
                    x.append(i-5)
    x.append(last-5)
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
    #delta=3
    #x = range(0,last,delta)
    #coordination_shell(VW, L_cont,x,sigma=1.)
    #coordination_shell(VW, L_cont,x,sigma=1.5)
    #coordination_shell(VW, L_cont,x,sigma=2.)
    #coordination_shell(VW, L_cont,x,sigma=2.5)
    #coordination_gyration(VW, L_cont, x)
    #cubic_order_animate(VW,Z,L_cont,x)
    #linker_gauss_test(M,VW, Z,L, x, rcut=1.0, step=5e4)
    #linker_gauss(M,VW, Z,L, x, rcut=1.0, step=5e4)
    #cubic_order(VW,Z,L_cont,x)
    #cube_connected_rotate(VW,Z,M, L_cont)
    #delta = 10
    #x = range(0,VW.shape[0],delta)
    #x = [last-10,last-1]
    #######################
    ######################
    #msd(VW,L_cont[0],step=1)
    import MD.analysis.cubatic_order as c
    #polymer_end_vs_lattice(M,VW, L_cont, x)
    #polymer_interpenetrate(M,VW, Z,L_cont,x)
    #polymer_azimuth(M,VW, Z,L_cont, x)
    #nn_distance(VW, L_cont, x, rcut=25.0)
    #x = [80,586,793,1008]
    #x = [20,350]
    #draw_polymer(M,L_cont, x,NP=[0,117,2])
    #polymer_angle(M,VW, L_cont, x)
    #t = last-10
    #dt=5
    #x = [0]
    #VW = Average(VW,L_cont,t,dt,save='star.xyz')
    #neighbor_no_gauss(np.array([VW]), L_cont, x,count=12 ,step=5e4,rcut=35)
    #x = range(0,last-5,5)
    #polymer_gauss(M,VW, Z,L_cont, x)
    #for i in x:
    #    y = range(i,i+10)
    #    polymer_inter_gauss(M,VW, L_cont, y)
    #ssDNA_gauss(M,VW, Z,L_cont, x)
    #c.gaussmap_cluster(VW,Z,L_cont,x,1,scale=12)
    ##t = 400
    #dt = 10
    #VW = Average(VW,L_cont,t,dt,save='star.xyz')
    #Z = Average(Z,L_cont,t,dt,save='Z.xyz')
    #c.gaussmap(np.array([VW]),np.array([Z]),L_cont,[0],1)
    #neighbor_no_gauss(np.array([VW]), L_cont, [0],count=12, step=5e4,rcut=25)
#################################################
def run_compress():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    delta = 1
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
        VW = util.pickle_load('VW.pkl')
    except:
         VW=M.cord_auto(['V','W'])
         util.pickle_dump(VW,'VW.pkl')
    x = [0]
    count = 1
    lcut = 0
    L_length = [[L_cont[0][0]]]
    L_last = L_cont[0]
    for i in range(1,len(L_cont),delta):
        if L_cont[i] != L_last :
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
            x[-1] += i[1]
            count+=i[1]
    print x
    print L_length
    msd_compress(VW,L_cont,x,step=1)
if __name__ == '__main__':
    run_single()
    #run_compress()

