## \package MD.analysis.bond_order
# \brief methods for local bond ordering of a group of particles in
# simulation.
#
#

import numpy as np
from MD.base.spherical import xyz_to_spherical,qlm,ql,Sij
from MD.analysis.nearest_neighbor import nearest_neighbors_index

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
def solid_particles(M,L,rcut=17.5,l=6,c_cut=0.15,crystal=6,verbose=False):
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
    S=[]
    for i in range(M.shape[0]):
        S.append(Sij(Qlm[i],Qlm[neighbors[i][0]]))
    #for each particle lets counts its neighbors
    #if it has more than ccut we say it is a solid
    iscrystal=[0.0 for i in range(M.shape[0])]
    counter=0
    store_j = []
    for i in S:
        #if verbose:
        #    print i
        count=0
        for j in i:
            #this is the cutoff for declaring 2 particles connected
            if j.real>c_cut:
                store_j.append(j.real)
                count+=1
        if count >= crystal:
            iscrystal[counter]=1
        counter+=1
    try:
        print sum(store_j)/len(store_j)
    except:
        print 'no bonds'
    num=iscrystal.count(1)
    print 'crystals',num
    return num, iscrystal
