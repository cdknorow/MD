# -*- coding: utf-8 -*-
import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
import MD.plot.pyplot_eps as pyplot_eps
from MD.analysis.particle_distance import particle_distance
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.particle_distance import min_point_distance
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import second_nearest_neighbors_index
from MD.analysis.nearest_neighbor import count_neighbors_index
from MD.analysis.nearest_neighbor import min_distance
from MD.plot.histogram import histogram
from MD.plot.histogram import histogram_normal
from MD.plot.histogram import histogram_reg
from drift import drift_remove
import pickle
# Run from inside a folder
#########
reload(MD)
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
###################################################
# Find Static Structure Facto
# S is sf peaks with the noise removed
# Q is the corresponding q-vectors magnitude that go with the sf peaks
# Primitive vectors is the correspoinding q-vector
#####################################################
def structure_factor(VW, L, n_frames=10, n_start=1, 
        filters=0.05, dh=0.05,save='sf',l=10):
    import MD.analysis.sfactor as sf
    if os.path.isfile('qlattice.pkl'):
        crystal = pickle.load(open('qlattice.pkl','r'))
        return crystal
    #Find Structure Factor
    A = Average(VW,L,n_start,n_frames)
    return A
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
    #pyplot.plot_sf(Q, S,0.552381, linestyle='', marker='x', 
    #        xlabel=xlabel, ylabel=ylabel, save=save,xlen=20)
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
        if raw_input('Reselect? (y)/N') == 'N':
            con = False
            crystal = mb.make_bcc(a1,a2,a3,A[0],L[n_start],S=30,name='qlattice0.xyz')
            pickle.dump(crystal[0],open('qlattice0.pkl','w'))
            crystal = mb.make_bcc(a1,a2,a3,A[1],L[n_start],S=30,name='qlattice1.xyz')
            pickle.dump(crystal[0],open('qlattice1.pkl','w'))
            crystal = mb.make_bcc(a1,a2,a3,A[2],L[n_start],S=30,name='qlattice2.xyz')
            pickle.dump(crystal[0],open('qlattice2.pkl','w'))
            crystal = mb.make_bcc(a1,a2,a3,A[3],L[n_start],S=30,name='qlattice.xyz')
            pickle.dump(crystal[0],open('qlattice.pkl','w'))
    return crystal[0]
#Find the gaussian distribution around the ideal lattice sites
def gaussian(VW,crystal,L,frames,step): #get nearest neighbor distance to point d_x = []
    print step
    d_y = []
    d_z = []
    d_x = []
    index = []
    delta = .01
    minbin = 50
    for i in VW[frames[0]]:
        index.append(min_point_distance(i,crystal,L)[2])
    for k in frames:
        for j,i in enumerate(VW[k]):
            dist = points.dist(i,crystal[index[j]],L)[1]
            d_x.append(dist[0])
            d_y.append(dist[1])
            d_z.append(dist[2])
    #cont = True
    #while cont == True:
        #delta = float(raw_input("delta = "))
        #minbin = float(raw_input("bin = "))
    bins = np.arange(-minbin,minbin+delta,delta)
    bin_x, hist_x, max_hist = histogram_normal(d_x,bins)
    bin_y, hist_y, max_hist = histogram_normal(d_y,bins)
    bin_z, hist_z, max_hist = histogram_normal(d_z,bins)
    pickle.dump([[bin_x,hist_x],[bin_y,hist_y],[bin_z,hist_z]],open('gaussian%i.pkl'%step,'w'))


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
    try:
        V = util.pickle_load('V.pkl')
        if V.shape[0] != last:
            print 'generating new pickle'
            fail
    except:
        #V=M.cord_auto(['V','W'])
        V=M.cord_auto(['V'])
        if V.shape[1] < 1:
            V=M.cord_auto(['A0','A1','A'])
        V = drift_remove(V,L)
        util.pickle_dump(V,'V.pkl')
    x = [-30]
    print 'getting lattice'
    for i in x:
        crystal = structure_factor(V, L_cont, n_start=i,save='sf'+str(i),
                n_frames = 20,l=20, filters=0.005)
    print 'finding gaussian'

    delta = int(raw_input('delta ='))
    r = last / delta
    for i in range(0,last-(r-1),r):
        x = range(i,i+r)
        gaussian(V,crystal,L,x,step=i)

###for single directories
if __name__ == '__main__':
        run_()

