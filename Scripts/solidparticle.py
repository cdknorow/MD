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

reload(av)

###########################################################
# get the steinhardt parameters for a partticle in a system system
#########################################################
def steinhardt(VW,L,x,avg=1,step=5e4,rcut=30):
    print "Particle Q2_avg  Q4_avg  Q6_avg Q8_avg Q10_avg"
    for k in x:
        A = av.average_position(VW,k,avg,L)
        Q2_avg = (bond_order.steinhardt_order(A,L[k],l=2,rcut=rcut))
        Q4_avg = (bond_order.steinhardt_order(A,L[k],l=4,rcut=rcut))
        Q6_avg = (bond_order.steinhardt_order(A,L[k],l=6,rcut=rcut))
        Q8_avg = (bond_order.steinhardt_order(A,L[k],l=8,rcut=rcut))
        Q10_avg = (bond_order.steinhardt_order(A,L[k],l=10,rcut=rcut))
        for i in range(len(Q2_avg)):
            print i, Q2_avg[i], Q4_avg[i], Q6_avg[i], Q8_avg[i], Q10_avg[i]



#########################################################
#Find the number of solid particles within the system
#########################################################
def solid_particles(VW, L, VW_names, x, avg=1, step=5e4, bcc = True, sc = False, GenLattice = False):
    solid_particles = []
    solid_crystals = []
    try:
        solid_crystals = util.pickle_load('spcrystal.pkl')
        solid_particles = util.pickle_load('solidp.pkl')
    except:
        for k in x:
            print 'step', k
            A = av.average_position(VW,k,avg,L)
            if bcc:
                num, crystals = (bond_order.solid_particles(A,
                                L[k],c_cut=0.15, count=8,l=6,crystal=7))
            if sc:
                num, crystals = (bond_order.solid_particles(VW[k],
                                L[k],c_cut=0.175, count=6,l=4,crystal=5))
            solid_crystals.append(crystals)
            solid_particles.append(num)
        util.pickle_dump(solid_crystals,'spcrystal.pkl')
        util.pickle_dump(solid_particles,'solidp.pkl')

    fid = open('solidpart.txt','w')
    for i in range(len(x)):
        fid.write('%i %i %.2f\n'%(x[i], solid_particles[i], L[x[i]][0]))
    fid.close()

    if GenLattice:
        # make an xyz file of the crystals and there locations
        f = open('solidpcrystals.xyz','w')
        count=0
        for k in x:
            f.write(('%i \n\n')%(VW.shape[1]))
            for i in range((len(solid_crystals[count]))):
                if solid_crystals[count][i] == 1:
                        f.write(('%c %.2f %.2f %.2f \n')%('V',VW[k][i][0],VW[k][i][1],VW[k][i][2]))
                else:
                        f.write(('%c %.2f %.2f %.2f\n')%('W',4*L[k][0],4*L[k][0],4*L[k][0]))
            count += 1
        f.close()
    return x, solid_particles



if __name__ == '__main__':
    import MD.unit.make_lattice as make_bcc
    a = 17.
    basis = np.array([0.0,0.0,0.0])
    a1 = np.array([-0.5*a,0.5*a,0.5*a])
    a2 = np.array([0.5*a,-0.5*a,0.5*a])
    a3 = np.array([0.5*a,0.5*a,-0.5*a])
    L = [np.array([a*5,a*5,a*5])]
    VW, VW_names = make_bcc.make_lattice(a1,a2,a3,basis,L[0],S=30,name='bcc_unit_large.xyz')
    reload(bond_order)
    steinhardt(np.array([VW]),L,[0],rcut=30)
    #solid_particles(np.array([VW]), L, VW_names, [0])
    #solid_particles(np.array([VW]), L, VW_names, [0], bcc=False, sc=True)
