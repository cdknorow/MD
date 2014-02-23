
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
#get connections
def connections(M):
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
############################################
## returns the arrays with the drift removed
############################################
def drift_remove(VW,L):
    from MD.analysis.drift_remove import eliminate_drift
    VW = eliminate_drift(VW,L)
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
def write_xyz(V):
    fid = open('V.xyz','w')
    for k in range(V.shape[0]):
        fid.write('%i\n\n'%V.shape[1])
        for i in range(V.shape[1]):
            fid.write('V %.3f %.3f %.3f \n'%(
                V[k][i][0], V[k][i][1], V[k][i][2]))
    fid.close()
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
    last = M.frames
    try:
        #V = util.pickle_load('V.pkl')
        LB = util.pickle_load('LB.pkl')
    except:
    #    #V=M.cord_auto(['V'])
         LB=M.cord_auto(['LB'])
         util.pickle_dump(LB,'LB.pkl')
    drift_remove(LB,L)
if __name__ == '__main__':
    run_single()

