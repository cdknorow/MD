#this script is used to find the diffusion of groups of particles within a
#simulation
import numpy as np
from MD.plot.histogram import histogram_reg
from MD.analysis.nearest_neighbor import nearest_neighbors_index
from MD.analysis.nearest_neighbor import close_neighbors_point
import MD.analysis.bond_order as bond_order
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.msd import msd_no_drift
import MD.analysis.average as av
import MD.plot.pyplot as pyplot
import MD.util as util
import MD.plot.pyplot_eps as pyplot_eps
reload(av)
#used to find the movement of the crystal to the surface
def crystal_movement(VW,surface,L):
    skip = VW.shape[0]/len(surface)
    msd = []
    for i in range(10,len(surface[0])):
        print i
        for k in range(len(surface)):
            if surface[k][i] == 1:
                count = 0
                if surface[k-1][i] == 0:
                    if k+count >= len(surface)-5:
                        still_surface = False
                    else:
                        still_surface = True
                    while still_surface == True:
                        if surface[k+count][i]:
                            count += 1
                        elif surface[k+count+1][i]:
                            count += 2
                        elif surface[k+count+2][i]:
                            count += 3
                        elif surface[k+count+3][i]:
                            count += 4
                        else:
                            still_surface = False
                        if k+count >= len(surface)-5:
                            still_surface = False
                        #find the displacment of the particle from k to count
                    start = (k-count)*skip
                    finish = k*skip
                    if finish - start <=1:
                        finish += 1
                    spheres = np.zeros((finish-start,1,3))
                    for f in range(finish-start):
                        spheres[f] = VW[f*skip+k*skip][i]
                    x, MSD = msd_no_drift(spheres,L)
                    msd.append([MSD,k])
                    k = k+count
    frame = np.zeros((len(surface),1))
    MSD = np.zeros((len(surface),1))
    for i in msd:
        for j in range(len(i[0])):
            frame[i[1]+j] += 1
            MSD[i[1]+j] += i[0][j]
    for i in range(MSD.shape[0]):
        MSD[i] = MSD[i] / frame[i]
    print MSD
    return MSD
#used to find the jump distance of particles from gel to surface
def crystal_jump(VW,surface,L,time,delta):
    #skip = VW.shape[0]/len(surface)
    skip = 1
    jump_t = []
    jump_d = []
    lfilter = 10
    for i in range(len(surface[0])):
        for k in range(time,time+delta):
            if surface[k][i] == 1:
                if surface[k-1][i] == 0:
                    count_before = 0
                    count_after = 0
                    for j in range(lfilter):
                        count_after += surface[k+j][i]
                    for j in range(lfilter):
                        count_before += surface[k-j][i]
                    if (count_before <=2) and (count_after >=8):
                        spheres = np.zeros((50,1,3))
                        for f in range(k-50,k):
                            spheres[f-k] = VW[f*skip][i]
                        x, MSD = msd_no_drift(spheres,L)
                        jump_d.append(MSD[-1])
                        break
    print 'number of counts',len(jump_d)
    return jump_d
#used to find the msd of the crystal
def crystal_msd(VW,surface,L):
    #skip = VW.shape[0]/len(surface)
    skip = 1
    index = []
    for i in range(len(surface[0])):
        count = 0
        for t in range(15):
            count += surface[t][i]
        if count >= 11:
            index.append(i)
    spheres = np.zeros((len(surface),len(index),3))
    for f in range(len(surface)):
        for i in range(len(index)):
            spheres[f][i] = VW[f*skip][index[i]]
    x, MSD = msd_no_drift(spheres,L)
    pyplot.plot(x,MSD,xlabel='time',ylabel='msd cn',save='msd_critical')
#find the msd of particles with index
def index_msd(VW,index,L,t,save='msd_index'):
    #skip = VW.shape[0]/len(surface)
    skip = 1
    spheres = np.zeros((t,len(index),3))
    for f in range(t):
        for i,j in enumerate(index):
            spheres[f][i] = VW[f*skip][j]
    x, MSD = msd_no_drift(spheres,L)
    pyplot_eps.plot(x,MSD,xlabel='t',ylabel=r'$\sqrt{<(r(t)-r(0))^2>}$',save=save)
    return x,MSDj
#another method to find the jump distance
def projection_jump(VW,surface,proj,L,k,save='proj_hist'):
    #skip = VW.shape[0]/len(surface)
    skip = 1
    jump_d = []
    print len(surface)
    print VW.shape
    print proj.shape
    try:
        for i in surface:
            n, dist = close_neighbors_point(proj,VW[i],L)
            jump_d.append(dist)
    except:
        for i in surface:
            n, dist = close_neighbors_point(proj,VW[i-VW.shape[0]],L)
            jump_d.append(dist)
    print jump_d
    hist_s,xs=histogram_reg(jump_d,bins=10)
    pyplot.plot(xs,hist_s,xlabel='distance',ylabel='count',save=save)
    return jump_d
