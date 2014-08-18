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
import h5py
############ 
# Run from inside a folder
#########
reload(MD)
reload(pyplot)


#\brief find the polymer strething distances
def polymer_end_2_end(M, L, frames ):
    hpy = h5py.File('polymer_stetching.hdf5','w')
    delta = 1

    P = M.cord_frames(['M'],frames)
    V = M.cord_frames(['V'],frames)
    polymers = P.shape[1]/(2*V.shape[1])

    for count, k in enumerate(frames):
        counter = 0
        ps2pe = np.zeros(V.shape[1]*polymers)
        c2pe = np.zeros(V.shape[1]*polymers)
        c2ps = np.zeros(V.shape[1]*polymers)
        for i in range(V.shape[1]):
            for j in range(0,2*polymers,2):
                ps2pe[counter] = points.dist(P[count][i*(2*polymers)+j],P[count][i*(2*polymers)+j+1],L[k])[0]
                c2ps[counter] = points.dist(V[count][i],P[count][i*(2*polymers)+j],L[k])[0]
                c2pe[counter] = points.dist(V[count][i],P[count][i*(2*polymers)+j+1],L[k])[0]
                counter+=1

        dset = hpy.create_dataset("%ips2pe"%k, data = ps2pe)
        dset.attrs.create("N_atoms",V.shape[1]*polymers)
        dset.attrs.create("L",L[k][0])
        dset.attrs.create("n_frames",delta)
        dset = hpy.create_dataset("%ic2ps"%k, data = c2ps)
        dset.attrs.create("N_atoms",V.shape[1]*polymers)
        dset.attrs.create("L",L[k][0])
        dset.attrs.create("n_frames",delta)
        dset = hpy.create_dataset("%ic2pe"%k, data = c2pe)
        dset.attrs.create("N_atoms",V.shape[1]*polymers)
        dset.attrs.create("L",L[k][0])
        dset.attrs.create("n_frames",delta)
    hpy.close()


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
    last = M.frames
    x, xcount = util.get_box_size_steps(L_cont)
    polymer_end_2_end(M, L_cont, [50,70,90,110,130])

#For multiple directories
#if __name__ == '__main__':
#    y = []
#    [y.append(x[0]) for x in os.walk(os.getcwd())]
#    del y[0]
#    for i in y:
#        print i.split('/')[-1]
#    for directory in y:
#        if os.path.isdir(directory):
#            os.chdir(directory)
#            try:
#                run_single()
#            except:
#                'directory failed',directory
#                pass
#            print '##\nfinished with:',directory

if __name__ == '__main__':
    run_single()
#    #run_compress()
