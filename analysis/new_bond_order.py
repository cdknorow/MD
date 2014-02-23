
## \package MD.analysis.bond_order
# \brief methods for local bond ordering of a group of particles in
# simulation.
#
#

import numpy as np
from MD.base.spherical import xyz_to_spherical,qlm,ql,Sij
from MD.analysis.nearest_neighbor import count_neighbors_index
from MD.analysis.nearest_neighbor import cn_tree_index
from MD.base.octree import node

## \brief Uses spherical harmonics to identify solid particles 
# a specfic frame.
#
# \return num number of solid particles in group of particles
#
# \param M a group of atoms in matrix form [atoms][x,y,z]
# \param L box-periodicty
# \param r_cut cutoff for neighbor search
# \param l bond order,6 for bcc, 4 for simple cubic 
# \param c_cut crystal value cutoff for identifying solid particle
# \param crystal number of surrounding crystal particles to call a particle a crystal 
#############################
def Qmap(M,L,count=14,l=6,c_cut=0.15,crystal=6,verbose=False,rcut=30):
    #Find nearest neighbors for all particles and convert to spherical
    #######
    #neighbors[i][0]=index of neighbors
    #neighbors[i][1]=array[i][x,y,z] of distance between neighbors
    ########
    neighbors=[]
    print 'finding octree'
    print L
    otree = node(L,delta=8)
    otree.add_particles(M)
    print 'octree created'
    for i in range(M.shape[0]):
        N = cn_tree_index(M,i,otree,L,count=count,rcut=rcut)
        neighbors.append([N[0],xyz_to_spherical(N[1])])
        #include itself
        neighbors[-1][0].append(i)
        if verbose:
            print 'neighbors',len(N[0])
    #find qlm of each particle 
    #######
    #store as values m=-l-l in nparray
    #Qlm[i]=np.array(ql-m:ql+m) for particle i
    #Qlm[i]=qlm(i)
    ######
    print M.shape
    Ql=np.zeros((M.shape[0]),dtype=complex)
    qlm_avg=np.zeros((2*l+1),dtype=complex)
    for i in range(M.shape[0]):
        print i
        for m in range(-l,l+1):
            for k in neighbors[i][0]:
                qlm_avg[m+l]+=qlm(neighbors[k][1],l,m)
            qlm_avg[m+l]/=len(neighbors[i][0])
        Ql[i]=ql(qlm_avg,l)

    print Ql
    return Ql
## \brief Calculates Steihardt order parameters  
#
# \return num steihard order parameter for particles
#
# \param M a group of atoms in matrix form [atoms][x,y,z]
# \param L box-periodicty
# \param r_cut cutoff for neighbor search
# \param l bond order,6 for bcc, 4 for simple cubic 
# \param c_cut crystal value cutoff for identifying solid particle
# \param crystal number of surrounding crystal particles to call a particle a crystal 
#############################
def steinhardt_order(M,L,l=6,rcut=17.5,verbose=False):
    print 'l',l
    #Find nearest neighbors for all particles and convert to spherical
    #######
    #neighbors[i][0]=index of neighbors
    #neighbors[i][1]=array[i][x,y,z] of distance between neighbors
    ########
    neighbors=[]
    for i in range(M.shape[0]):
        N=nearest_neighbors_index(M,i,L,rcut=rcut)
        neighbors.append([N[0],xyz_to_spherical(N[1])])
        if verbose:
            print 'neighbors',len(N[0])
    #find qlm of each particle 
    #######
    #store as values m=-l-l in nparray
    #Qlm[i]=np.array(ql-m:ql+m) for particle i
    #Qlm[i]=qlm(i)
    ######
    Qlm=np.zeros((M.shape[0],2*l+1),dtype=complex)
    qi=np.zeros((1,2*l+1),dtype=complex)
    for i in range(M.shape[0]):
        for m in range(-l,l+1):
            qi[0][m+l]=qlm(neighbors[i][1],l,m)
        Qlm[i][:]=qi
    #for nearest neighbors find the Sij vector
    Ql=[]
    for i in range(M.shape[0]):
        Ql.append(ql(Qlm[i],l))
    print sum(Ql)/len(Ql)
    return Ql
