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
from MD.analysis.nearest_neighbor import count_neighbors_index
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
    out = open('msd.txt','w')
    for i in range(len(msd)):
        out.write('%i %.2f\n'%(x[i],msd[i]))
    out.close()
    pyplot.plot(x,msd,xlabel='Time',ylabel='msd',save=save)
    #pyplot.plot3(x,msd,x,(6*D*x)**0.5,x,(6*D2*x)**0.5,xlabel='Time',ylabel='msd',label1='msd',label2='%.2f'%D,label3='%.2f'%D2,save=save,showleg=True)
    #pyplot.plot3(x,msd,x,(6*D2*x)**0.5,x,msd2,xlabel='Time',ylabel='msd',label1='msd',label2='%.2f'%D2,label3='%.2f'%D,save=save,showleg=True)
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
## 	type ordering
############################################
def type_ordering(VW,L):
    N = VW.shape[1]*8
    order = []
    count = 0
    for k in range(VW.shape[0]):
        print k
        print count
        count = 0
        for i in range(VW.shape[1]):
            NN=count_neighbors_index(VW[k],i,L[k],count=8)[0]
            NN.sort()
            if i < VW.shape[1]/2:
                for j in NN:
                    if j >= VW.shape[1]/2:
                        count+=1
            if i >= VW.shape[1]/2:
                for j in NN:
                    if j < VW.shape[1]/2:
                        count+=1
        order.append(count/float(N))
    x = range(VW.shape[0])
    fid = open('typeorder.txt','w')
    for i in range(len(order)):
        fid.write(('%i %.2f\n'%(i,order[i])))
    fid.close()
    pyplot.plot(x,order,
            xlabel='time',ylabel='ordering',save='typeorder')

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
    Q, S, primitive_vectors = sf.sf_filter(Q, S, primitive_vectors,
            filters=filters, save = 'sf%i'%n_start)
    return Q, S, primitive_vectors
def run_binary():
    #print out directory dirname = os.getcwd().partition('/')[-1]
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    print L
    last = M.frames
    L_cont = M.box_volume()
    try:
        V = util.pickle_load('V.pkl')
        W = util.pickle_load('W.pkl')
        VW = util.pickle_load('VW.pkl')
    except:
        #delta = 25
        #V=M.cord_range(['V'],delta=delta)
        #W=M.cord_range(['W'],delta=delta)
        #VW=M.cord_range(['V','W'],delta=delta)
        V=M.cord_auto(['V'])
        W=M.cord_auto(['W'])
        VW=M.cord_auto(['V','W'])
        VW_index=M.get_names(['V','W'])
        #VW,V,W = drift_remove_all(VW,V,W,L,VW_index)
        util.pickle_dump(V,'V.pkl')
        util.pickle_dump(W,'W.pkl')
        util.pickle_dump(VW,'VW.pkl')
    delta = 30
    x = range(0,V.shape[0],delta)
    print V.shape[0]
    print x
    #plt.close()
    msd(VW,L,time_scale=delta)
    type_ordering(VW,L_cont)
    #plt.close
    #for i in x:
    #    print "finding s(q)"
    #    structure_factor(VW, L_cont[i], n_start=i,save='sf'+str(i))
    #    #structure_factor(V, L, n_start=i,save='sfV'+str(i))
    #    #structure_factor(W, L, n_start=i,save='sfW'+str(i))
    #    print "finding g(s)"
    #    plt.close
        #distance_distribution(V,W,L_cont[i],n_start=i)
    #    plt.close()
    #    #plt.close()
#For multiple directories
#if __name__ == '__main__':
#    for f in sorted(os.listdir("./")):
#        if os.path.isdir(f):
#            os.chdir(f)
#            try:
#                run_binary()
#                print '##\nfinished with:',f
#            except:
#                print 'failed in this directory'
#            os.chdir('../')
if __name__ == '__main__':
##    #run_debug()
##    #run_all()
##    #run_single()
##    #run_simple()
    run_binary()

