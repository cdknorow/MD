## \package MD.analysis.graph 
# \brief This module is used to find a network graph between particle over a series
# of frames using the package python-graph
# 
# dependencies
# python-graph
#
# This code is a method for finding the number of clusters of networks during a
# simulation where Major particles with many connecting sites are able to form
# networks with other major particles through the connecting sites. 
# Here we analyze the number of clusters formed through these networks over time.

import networkx as nx

## \brief find the degree of connectivity among the graph
#
#
# \param gr graph
# 
# \returns count where each index represents the number of edges 
#   and each index contains the number of atoms with that amount of edges
#
def degree_of_connection(gr):
    T = len(gr.adjacency_list())
    count = [0,0,0,0,0,0,0,0,0,0]
    for i in range(T):
        for j in gr.neighbors(i):
            try:
                count[gr.number_of_edges(i,j)] += 1
            except:
                pass
    return count

## \brief find the total number of neighbors for each and sum over all
#
#
# \param gr graph
# 
# \returns total number of neighbors in gr
#
def neighbor_find(gr):
    count = 0
    for i in gr.nodes():
        count += len(gr.neighbors(i))
    return count
## \brief create a graph using networkx package 
#
# \returns number of independent sets at each frame
# \returns list of connected points in each frame
#
# \param connections connecting points at each frame
# \param N_nodes the number of central nodes
# \param vertex_points the number of vertex points at each node
def grapher(connections, N_nodes, vertex_points = 1):
    networks = []
    degree = []
    neighbors = []
    num_neighbors = []
    for i in connections:
        gr = nx.MultiGraph()
        gr.add_nodes_from(range(0,N_nodes))
        for j in i:
            #find the two connections
            #we need to divide by the number of sites to get the right
            #major particle index
            A=j[0]/vertex_points
            B=j[1]/vertex_points
            #add the connections to a list
            gr.add_edge(A,B)
        networks.append(nx.connected_components(gr))
        degree.append(degree_of_connection(gr))
        num_neighbors.append(neighbor_find(gr))
        neighbors.append(gr.adjacency_list())
    num_networks = []
    deg = [ [] for x in range(len(degree[0])) ]
    for i in networks:
        num_networks.append(len(i))
        deg[0].append(len(i))
    for i in degree:
        for j in range(1,len(i)):
            deg[j].append(i[j]/2)
    return networks, num_networks, deg, neighbors, num_neighbors, gr

## \brief create a graph using networkx package 
#
# \returns number of independent sets at each frame
# \returns list of connected points in each frame
#
# \param connections connecting points at each frame
# \param N_nodes the number of central nodes
# \param vertex_points the number of vertex points at each node
def grapher_2(conA, conB, N_nodes = 54, N_linker = 500, vertex_points = 25):
    networks=[]
    nA = vertex_points * N_nodes/2
    nB = vertex_points * N_nodes/2
    connections=[]
    for i in range(len(conA)):
        gr = nx.MultiGraph()
        gr.add_nodes_from(range(0,nA+nB+N_linker))
        for j in conA[i]:
            A=j[0]
            B=j[1]+nB
            #add the connections to a list
            gr.add_edge(A,B)
        for j in conB[i]:
            A=j[0]+nA
            B=j[1]+nA
            #add the connections to a list
            gr.add_edge(A,B)
        #This will give a graph of connections
        con = nx.connected_components(gr)
        #remove the linker component from the graph
        for n in range(len(con)):
            for p in range(nA+nB,nA+nB+N_linker):
                if con[n].count(p)==1:
                    con[n].remove(p)
        #filter empty sets
        list2 = [x for x in con if x != []]
        #filter sets with only 1 value
        list3 = [x for x in list2 if len(x) != 1]
        connections.append(list3)

    #Now find the graph between A and B particle
    networks, num_networks = grapher(connections, N_nodes, vertex_points =
            vertex_points)
    return networks, num_networks

def grapher_2(conA, N_nodes = 54, N_linker = 500, vertex_points = 25):
    networks=[]
    nA = vertex_points * N_nodes
    connections=[]
    for i in range(len(conA)):
        gr = nx.MultiGraph()
        gr.add_nodes_from(range(0,nA+nB+N_linker))
        for j in conA[i]:
            A=j[0]
            B=j[1]+nA
            #add the connections to a list
            gr.add_edge(A,B)
        #This will give a graph of connections
        con = nx.connected_components(gr)
        #remove the linker component from the graph
        for n in range(len(con)):
            for p in range(nA+nB,nA+nB+N_linker):
                if con[n].count(p)==1:
                    con[n].remove(p)
        #filter empty sets
        list2 = [x for x in con if x != []]
        #filter sets with only 1 value
        list3 = [x for x in list2 if len(x) != 1]
        connections.append(list3)

    #Now find the graph between A and B particle
    networks, num_networks = grapher(connections, N_nodes, vertex_points =
            vertex_points)
    return networks, num_networks


## \brief create a second neighbor graph using networkx package 
#
# \returns number of independent sets at each frame
# \returns list of connected points in each frame
#
# \param connections connecting points at each frame
# \param N_nodes the number of central nodes
# \param vertex_points the number of vertex points at each node
def second_grapher(connections, N_nodes, vertex_points = 1,
        f1cut = 1, f2cut =1, greater = False):
    networks = []
    degree = []
    neighbors = []
    num_neighbors = []
    gr = nx.MultiGraph()
    gr.add_nodes_from(range(0,N_nodes))
    for i in connections:
        #find the two connections
        #we need to divide by the number of sites to get the right
        #major particle index
        A=i[0]/vertex_points
        B=i[1]/vertex_points
        #add the connections to a list
        gr.add_edge(A,B)
    gsmall = nx.MultiGraph()
    degree.append(degree_of_connection(gr))
    #filters out a unique set in array
    def filter(seq):
        # not order preserving 
        set = {}
        map(set.__setitem__, seq, [])
        return set.keys()
    sn = []
    #find a list of second order neighbor connectsion 
    #with a defree of connections greater than f1/f2cut
    if greater:
        for i in range(N_nodes):
            nn = gr.neighbors(i)
            second_neighbors = []
            for j in nn:
                if gr.number_of_edges(i,j) >= f1cut:
                    for k in gr.neighbors(j):
                        if gr.number_of_edges(j,k) >= f2cut and k != i:
                            second_neighbors.append(k)
            sn.append(filter(second_neighbors))
    else:
        for i in range(N_nodes):
            nn = gr.neighbors(i)
            second_neighbors = []
            for j in nn:
                if gr.number_of_edges(i,j) == f1cut:
                    for k in gr.neighbors(j):
                        if gr.number_of_edges(j,k) >= f2cut and k != i:
                            second_neighbors.append(k)
            sn.append(filter(second_neighbors))
    return sn

## \brief create a graph using networkx package 
#
# \returns number of faces between two cubes which have connections
#
# \param connections connecting points at each frame
# \param N_nodes the number of central nodes
# \param vertex_points the number of vertex points at each node
def sides(connections, N_nodes, vertex_points = 1, cut = 5):
    #Cube = np.zeros((len(connections),N_nodes,N_nodes,6))
    Cube = [[[[0.0 for i in range(6)] for j in range(N_nodes)] for
        k in range(N_nodes)] for l in range(len(connections))]
    for k,i in enumerate(connections):
        for j in i:
            #find the two connections
            #we need to divide by the number of sites to get the right
            #major particle index
            A=j[0]/vertex_points
            B=j[1]/vertex_points
            side1 = j[0]%vertex_points/(vertex_points/6)
            side2 = j[1]%vertex_points/(vertex_points/6)
            Cube[k][A][B][side1]+=1
            Cube[k][B][A][side2]+=1
            #add the connections to a list
    count = 0
    counter = []
    for i in range(len(Cube)):
        count = 0
        for j in range(len(Cube[i])):
            for k in range(len(Cube[i][j])):
                if Cube[i][j][k].count(0) == cut :
                    count += 1
        counter.append(count/2)

    print counter
    return counter

## \brief create a graph using networkx package 
#
# \returns number of faces between two cubes which have connections
#
# \param connections connecting points at each frame
# \param N_nodes the number of central nodes
# \param vertex_points the number of vertex points at each node
def sides_mod(connections, N_nodes, vertex_points = 1, cut = 5):
    #Cube = np.zeros((len(connections),N_nodes,N_nodes,6))
    Cube = [[[[] for i in range(6)] for j in range(N_nodes)] 
        for l in range(len(connections))]
    for k,i in enumerate(connections):
        for j in i:
            #find the two connections
            #we need to divide by the number of sites to get the right
            #major particle index
            A=j[0]/vertex_points
            B=j[1]/vertex_points
            side1 = j[0]%vertex_points/(vertex_points/6)
            side2 = j[1]%vertex_points/(vertex_points/6)
            if Cube[k][A][side1].count(B) == 0:
                Cube[k][A][side1].append(B)
            if Cube[k][B][side2].count(A) == 0:
                Cube[k][B][side2].append(A)
            #add the connections to a list
    count = 0
    counter = []
    for i in range(len(Cube)):
        count = 0
        for j in range(len(Cube[i])):
            for k in range(len(Cube[i][j])):
                if len(Cube[i][j][k]) == cut :
                    count += 1
        counter.append(count)

    print counter
    return counter
