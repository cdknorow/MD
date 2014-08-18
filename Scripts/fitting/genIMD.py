import sys
import os
import numpy as np
sys.path.append('/home/cdknorow/Dropbox/Software/')
import MD
from MD.analysis.nearest_neighbor import nearest_neighbors_index
import MD.analysis.bond_order as bond_order
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.msd import msd_no_drift
import MD.analysis.average as av
import MD.plot.pyplot as pyplot
import MD.util as util
import MD.base.points as points

def gen_config(delta=5):
    #setup
    dirname = os.getcwd().partition('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    last = M.frames
    V=M.cord_auto(['V','W'])
    if V.shape[1] < 1:
        V=M.cord_auto(['A'])
    x = range(0,last-delta,delta)
    for i in x:
        A = Average(V,L_cont,i,delta)
        util.print_xyz(A,L[0],atoms=V.shape[1],save='config.%i.dat'%i)
    fid = open('center.xyz','w')
    fid.write('1\n\n')
    fid.write(('V %.2f %.2f %.2f'%(L[0]/2.,L[1]/2.,L[2]/2.)))
    fid.close()
