## \package MD.analysis.percolate 
# \brief This module is used to find a network graph between particle over a series
# of frames
# 
# This code is a method for finding the number of clusters of networks during a
# simulation where Major particles with many connecting sites are able to form
# networks with other major particles through the connecting sites. 
# Here we analyze the number of clusters formed through these networks over time.

from MD.base.points import unique_array, intersection, unique, union

## \internal
# \brief search function for percolated
#
# \param connected array of connections see connections
# \param A index of particle A
# \param B index of particle B
def search(connected,A,B):
    add=False
    for i in range(len(connected)):
        #is it in there already
        if connected[i].count(A)==1:
            if connected[i].count(B)==0:
                connected[i].append(B)
        if connected[i].count(B)==1:
            if connected[i].count(A)==0:
                connected[i].append(A)
    #look through connectd if there is an intersection, take the union and
    #delete what was there
    for i in range(len(connected)):
        for j in range(len(connected)):
            if intersection(connected[i],connected[j])!=[]:
                connected[i]=union(connected[i],connected[j])
                connected[j]=connected[i]
    #Removed any doubles
    connected=unique(connected)
    return connected

## \brief find network graph between two types of particles which connect
# 
# \returns array of arrays containg paritcle index in each network
# [frame][[M1,M6,M8],[M2,M7]]...] where M1 would be the index of a major
# particle 
# \returns number of clusters in each frame
#
# \param network same as connections [frame][[p1,p6][p2,p9]] where p1 and p6 are
# connected etc. p1 refers to the index of the particle
# \param N_sites number of sites attatched to each major particle
# \param N_atoms number of major particles which can be connected too
#######################################################
def percolate(network,N_sites,N_atoms=54):
    #search for
    print "Finding networks"
    #### look through the network
    print len(network)
    network = unique_array(network)
    print len(network)
    print network
    percolated=[]
    for i in network:
        connected=[]
        #we need to make the initial have all separate connections
        #with each one being in a separate array
        for k in range(N_atoms):
            connected.append([k])
        for j in i:
            #find the two connections
            #we need to divide by the number of sites to get the right
            #major particle index
            A=j[0]/N_sites
            B=j[1]/N_sites
            #add the connections to a list
            connected=search(connected,A,B)
        percolated.append(connected)
    networks=[]
    for i in percolated:
        networks.append(len(i))
    print percolated
    #percolated is the data, networks is the number of networks in each frame
    return percolated,networks
