
import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
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
from MD.analysis.rotation import single_rotation_box









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
    #Average_position
    stmp,qx = sf.sfactor(np.array([VW]),L=L,l=l)
    S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
    #Plot the graph
    xlabel = '$|\\vec{q}$ $|$'
    ylabel = '$S(\\vec{q}$ $)$ '
    pyplot.plot(Q, S, linestyle='', marker='x', 
            xlabel=xlabel, ylabel=ylabel, save=save)
    util.pickle_dump(Q,'Vsfx%i.pkl'%n_start)
    util.pickle_dump(S,'Vsfy%i.pkl'%n_start)
    print 'finished reconstruction'
    Q, S, primitive_vectors = sf.sf_filter(Q, S, primitive_vectors,
            filters=filters,save=save)
    print 'sorting peaks'
    return Q, S, primitive_vectors


fid = open('qlattice.xyz','r')
N_atoms = int(fid.readline())
L = float(fid.readline().split('=')[-1])
L_cont = [L,L,L]
M = fid.readlines()
V = []
count = 0
index = []
for line in M:
    index.append(line.split()[0])
    V.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
fid.close()
structure_factor(V, L_cont, n_start=0,save='sf_qlattice',n_frames = 1,l=20)
