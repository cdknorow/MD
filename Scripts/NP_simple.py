
# -*- coding: utf-8 -*-
import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import time
import numpy as np
import math
import matplotlib.pyplot as plt
import h5py
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
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
    print x
    print msd
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
def Average(VW,L,n_start,n_frames):
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
    fid = open('star.xyz','w')
    fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
    for i in range(Average.shape[0]):
        fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
    fid.close()
    return Average
############################################
## returns the average position over a few frames
############################################
def average_run(VW,L,delta):
    #fAverage_position
    Average = np.zeros((VW.shape[0]/delta,VW.shape[1],VW.shape[2]))
    for k in range(VW.shape[0]/delta):
        base=k*delta
        for d in range(delta):
            step = d + base
            for i in range(VW.shape[1]):
                #check the distance between
                for j in range(3):
                    if abs(VW[step][i][j]-VW[base][i][j]) < L[step][j]/2:
                        Average[k][i][j] += VW[step][i][j]
                    elif VW[step][i][j]-VW[base][i][j] < -L[step][j]/2:
                        Average[k][i][j] += VW[step][i][j]+L[step][j]
                    elif VW[step][i][j]-VW[base][i][j] > L[step][j]/2:
                        Average[k][i][j] += VW[step][i][j]-L[step][j]
   # fix any points that may be outside of the box after averaging
    for k in range(Average.shape[0]):
        for i in range(Average.shape[1]):
            for j in range(3):
                Average[k][i][j] /= delta
                if Average[k][i][j] > L[k*delta][j]:
                    Average[k][i][j] -= L[k*delta][j]
                if Average[k][i][j] < -L[k*delta][j]:
                    Average[k][i][j] += L[k*delta][j]
    fid = open('star.xyz','w')
    for k in range(Average.shape[0]):
        fid.write('%i\nL=%.4f\n'%(Average.shape[1],L[k*delta][0]))
        for i in range(Average.shape[1]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[k][i][0],Average[k][i][1],Average[k][i][2])))
    fid.close()
    return Average
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
    #pyplot.plot(x, y, xlabel='Volume', ylabel=label, save=save)
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
    pyplot.plot(BB[0][1:],BB[1][1:],xlabel='s',ylabel='g(s)',save='nnplot_'+save_name)
    out = open('dist%i.txt'%n_start,'w')
    out.write('\n')
    for i in range(len(BB[1])):
        out.write(('%.5f %.5f\n')%(BB[0][i],BB[1][i]))
    out.close()
############################################
## Generate Pickle of Distances
############################################
def PMF(hpy, V,L,n_frames=15,n_start=1,end=1):
    start=n_start
    finish=n_start+n_frames
    time_start = time.clock()
    distance=particle_distance(V[start:finish],V[start:finish],L)
    print "g(r) Elapsed Time -" , time.clock()-time_start
    time_start = time.clock()
    dset = hpy.create_dataset("%i"%n_start, data = np.array(distance))
    dset.attrs.create("N_atoms",V.shape[1])
    dset.attrs.create("V",L[0]*L[1]*L[2])
    dset.attrs.create("L",L[0])
    dset.attrs.create("n_frames",n_frames)

    #util.marshal_dump([V.shape[1],L[0],n_frames,distance], 'pmf_dist%i.pkl'%n_start)
    print "Marshal Elapsed Time -" , time.clock()-time_start

#####################################################
# Find Static Structure Facto
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor(VW, L, n_frames=10, n_start=1, 
        filters=0.05, dh=0.05,save='sf',l=10):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    A = Average(VW,L,n_start,n_frames)
    #A = VW[n_start]
    stmp,qx = sf.sfactor(np.array([A]),L=L[n_start],l=l)
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
    fid = open('sftotal.txt','w')
    fid.write('S  Q\n')
    for i in range(len(S)):
        fid.write('%.5f %.5f\n'%(S[i],Q[i]))
    fid.close()
    #Plot the graph
    xlabel = '$|\\vec{q}$ $|$'
    ylabel = '$S(\\vec{q}$ $)$ '
    #sp4 = 0.634698
    #sp12 = 0.530493
    #sc 1.121997
    #pyplot_eps.plot_sf(Q, S,0.552381, linestyle='', marker='x', 
            #xlabel=xlabel, ylabel=ylabel, save=save,xlen=20)
    pyplot.plot_sf(Q, S,0.552381, linestyle='', marker='x', 
            xlabel=xlabel, ylabel=ylabel, save=save,xlen=20)
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
            filters=filters,save=save)
    import MD.analysis.reconstruct_lattice as rl
    import MD.unit.make_bcc as mb
    con = True
    while con:
        a1,a2,a3 = rl.reconstruct_lattice(Q,primitive_vectors,save+'recipricol.txt')
        print a1,a2,a3
        mb.make_bcc(a1,a2,a3,A[0],L[n_start])
        if raw_input('Reselect? (y)/N') == 'N':
            con = False
    return Q, S, primitive_vectors
#####################################################
# Find Self Intermediate Scattering Function
#
##S is sf peaks with the noise removed
##Q is the corresponding q-vectors magnitude that go with the sf peaks
##Primitive vectors is the correspoinding q-vector
#####################################################
def self_intermediate_scattering(VW, L, x, n_start=1, 
        filters=0.05, dh=0.05,save='sf',l=5,delta=5,calc=True):
    import MD.analysis.sfactor as sf
    #Find Structure Factor
    if calc:
        import time
        from MD.canalysis.int_scattering import int_scatter
        F = []
        for i in range(delta):
            print i
            print "caculating Average" 
            print VW.shape
            print "Structure Factor" 
            A = VW[n_start+i]
            stmp,qx = sf.sfactor(np.array([A]),L=L[n_start],l=l)
            S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
            S,Q,primitive_vectors = sf.sf_filter_max(S,Q,primitive_vectors,filters=0.5)
            index = S.index(max(S))
            k = primitive_vectors[S.index(max(S))]
            #take the average over different starting configurations
            print "Calculating Int Scattering" 
            start = time.clock()
            scat_c = int_scatter(VW[n_start+i:],np.array(L[0]),np.array(k),
                    np.array(x[:-(n_start+i)]))
            end = time.clock()
            print 'cython runtime'
            print end-start
            F.append(scat_c)
        F_avg = np.zeros((len(F[-1])))
        for i in F:
            for j in range(len(F[-1])):
                F_avg[j] += i[j]
        F_avg /= delta
        #reset the zero of x so that the plot starts at the correct time
        x = range(len(F[-1]))
        util.pickle_dump([x,F_avg],'si_scat.pkl')
    xlabel = 't'
    ylabel = 'F(q,t)'
    save = 'int_scatter'
    A = util.pickle_load('si_scat.pkl')
    x = A[0]
    F_avg = A[1]
    #Kohlraush-Williams-Watts function
    # A exp[-(t/tau_b)^beta]
    beta = 1.4
    tau_a = 6
    K = []
    K_x = []
    fid = open('self_scat.dat','w')
    for i in range(len(A[0])):
        fid.write('%i %.4f\n'%(A[0][i],A[1][i]))
    fid.close()
    for i in np.arange(1,len(x),0.1):
        K.append(F_avg[1] * math.exp(-(i/tau_a)**beta))
        K_x.append(i)
    #pyplot.plot(x,F_avg, linestyle='',marker='x', xlabel=xlabel,
    #        ylabel=ylabel,logx=True,limx=(0,100), save=save)
    #pyplot.plot2(x,F_avg,K_x,K, xlabel=xlabel,
            #ylabel=ylabel,logx=True,limy=(0,1.0),limx=(0,100), save=save)
######################
# FInd the solid particles in a simulation
#####################
def solid_particles_bcc(VW,L,x,step=5e4,rcut=False):
    import MD.dna_scripts.solidparticle as sp
    reload(sp)
    if os.path.exists('bccrystals.xyz') == False:
        sp.solid_particles_VW(VW,L,['V','V'],x,step=5e4,avg=2)
    #sp.crystal_size(L,VW.shape[0]/skip,VW.shape[1])
    #if os.path.exists('sccrystals.xyz') == False:
    #    sp.solid_particles_VW(VW,L,['V','V'],skip=skip,step=5e4,rcut=False)
    #sp.crystal_size_cube(L,VW.shape[0]/skip,VW.shape[1])

######################
# FInd the solid particles in a simulation
#####################
def solid_particles_sc(VW,L,x,step=5e4,rcut=15,delta=25):
    import MD.dna_scripts.solidparticle as sp
    reload(sp)
    if os.path.exists('sccrystals.xyz') == False:
        sp.solid_particles_sc(VW,L,['V','W'],x,step=5e4)
    sp.sc_surface(VW,L,rcut=rcut,delta=delta)
    #sp.crystal_size(L,VW.shape[0]/skip,VW.shape[1])
    #if os.path.exists('sccrystals.xyz') == False:
    #    sp.solid_particles_VW(VW,L,['V','V'],skip=skip,step=5e4,rcut=False)
    #sp.crystal_size_cube(L,VW.shape[0]/skip,VW.shape[1])



##\brief 
#################################################
# What to run
#################################################
def run_():
    #setp
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
    new = False
    try:
        V = util.pickle_load('V.pkl')
        if V.shape[0] != last:
            print 'generating new pickle'
            asdf
        #Z = util.pickle_load('Z.pkl')
    except:
        #V=M.cord_auto(['V','W'])
        new = True
        V=M.cord_auto(['V'])
        if V.shape[1] < 1:
            V=M.cord_auto(['A','A0','A1'])
        #Z=M.cord_auto(['Z'])
        #VW_index=M.get_names(['V','W'])
        #Z_index=M.get_names(['Z'])
        #V = drift_remove(V,L,VW_index)
        #Z = drift_remove(Z,L,Z_index)
        #dump_xyz(V)
        #dump_xyz(Z)
        util.pickle_dump(V,'V.pkl')
        #util.pickle_dump(Z,'Z.pkl')
    # if (os.path.isfile('msd.pkl') is False):
    #    msd(V,L)
    print V
    asdf
    x, xcount = util.get_box_size_steps(L_cont)
    print x
    x_set = []
    n_frames = []
    print xcount
    print x
    #F = filter(lambda x: x > 10, xcount)
    for i in range(len(xcount)):
        if xcount[i] > 10:
            x_set.append(x[i])
            n_frames.append(xcount[i]-10)
    print n_frames
    print x_set

    f = h5py.File('../pmf_dist.hdf5','w')
    for j,i in enumerate(x_set):
        if  i + n_frames[j] <= V.shape[0]:
            if L_cont[i] == L_cont[i+n_frames[j]]:
                print "L = ",L_cont[i][0]
                print "rho = ",V.shape[1]/L_cont[i][0]**3
                print "frame = ", i 
                print "n_frame", n_frames[j]
                PMF(f, V,L_cont[i],n_start=i,n_frames=n_frames[j])
            else:
                print "L is not the same reduce the n_frames"
    f.close()
    #x = [20]
    #end_end(M,L_cont)
    #cubics(M,VW,L,var['ndna'])
    #rotations(V,Z,L_cont,'Z')
    #plt.close()

    #plt.close()
    #delta = 25
    #x = range(0,V.shape[0],delta)
    #V = average_run(V,L_cont,10)
    #x = range(0,195,5)
    #x = [5]
    #solid_particles_bcc(V,L_cont,x)
    #Q6Q4_map(V,L_cont,x,rcut=15.5*2**0.5)
    #solid_particles_sc(V,L_cont,x)
    ##cube_angle(Z,V,L_cont)
    #x = range(0,V.shape[0],delta)
    #x = range(4500,4550)
    #x=[910]
    #x = [500,630]
    #distance_distribution(V,L_cont,x)
    #x = [4500]
    #print last
    #x = range(0,last-2,1)
    ##x = range(1,200,1)
    #if (os.path.isfile('si_scat.pkl') is False) or new:
    #    self_intermediate_scattering(V, L_cont, x,
    #            n_start=1,calc=True,delta=5)
    #else:
    #    print 'skiping self_intermediate'
    #V = Average(V,L_cont,185,5)
    #x= [last-10]
    #for i in x:
    #    structure_factor(V, L_cont, n_start=i,save='sf'+str(i),
    #            n_frames = 1,l=20, filters=0.005)
    #    #structure_factor(V, L_cont, n_start=i,save='sf'+str(i),
    #    #        n_frames = 5,l=15, filters=0.05)
    #    plt.close()
    #    print "finding g(s)"
         #distance_distribution(V,L_cont[i],n_start=i,n_frames=5)
    #    plt.close()
#For multiple directories
#if __name__ == '__main__':
#    y = []
#    [y.append(x[0]) for x in os.walk(os.getcwd())]
#    del y[0]
#    for i in y:
#        print i.split('/')[-1]
#    for directory in y:
#        if directory[-1] == 'V':
#            print directory
#            if os.path.isdir(directory):
#                os.chdir(directory)
#                try:
#                    run_()
#                except:
#                    'directory failed',directory
#                    pass
#                print '##\nfinished with:',directory
#For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir(os.getcwd())):
#        if os.path.isdir(f):
#            print f
#            os.chdir(f)
#            os.chdir('V')
#            run_()
#            print '##\nfinished with:',f
#            os.chdir('../../')
###for single directories
if __name__ == '__main__':
        run_()

