import numpy as np
from MD.analysis.nearest_neighbor import nearest_neighbors_index
import MD.analysis.bond_order as bond_order
from MD.analysis.particle_distance import min_particle_distance
from MD.analysis.msd import msd_no_drift
import MD.analysis.average as av
import MD.plot.pyplot as pyplot
import MD.util as util
reload(av)

###########################################################
# solid particles specifiallcy for bcc systems
#########################################################

#Find the number of crystals, the size of the three largest crystals, and the width of the
#largest crystal in the x, y, and z direction
def crystal_size(L,num_frames,num_particles,rcut=15):
    #load all crystal positions
    num_frames = 777
    #rigid
    rcut = 25
    fid = open('bcccrystals.xyz','r')
    M = np.zeros((num_frames,num_particles,4))
    for frame in range(num_frames):
        fid.readline()
        fid.readline()
        for i in range(num_particles):
            row = fid.readline().split()
            M[frame][i][0]=float(row[1])
            M[frame][i][1]=float(row[2])
            M[frame][i][2]=float(row[3])
            if row[0] == 'V':
                M[frame][i][3]=0
            if row[0] == 'W':
                M[frame][i][3]=1
    #Find conected crystals 
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    import networkx as nx
    try:
        networks = util.pickle_load('cnetwork.pkl')
    except:
        connections = con.connections_min_max(M,M,L,rcut=rcut,rmin=4)
        for i in range(len(connections)):
            for j in range(len(connections[i])):
                connections[i][j][1] -= num_particles
        #Find groups of crystals via percolation, ie use distance as connection
        networks, num_networks, deg, neighbors, num_n, gr =   graph.grapher(connections,num_particles,1)
        util.pickle_dump(networks,'cnetwork.pkl')
    #Once groups of crystals are identified, find the size of each crystal
    #find the largest crystal in each frame
    num_crystal = []
    max_crystal = []
    second_crystal = []
    third_crystal = []
    crystal_cords_frame = []
    for k,i in enumerate(networks):
        crystal = []
        crystaltwo = 0
        crystalthree = 0
        max_c = 1
        max_twoc = 1
        max_threec = 1
        count = 0
        for j in i:
            if len(j) > max_c:
                crystal = j
                max_c = len(j)
                crystal_cords = j
            if len(j) > 1:
                count += 1
        for j in i:
            if len(j) < max_c:
                if len(j) > max_twoc:
                    crystaltwo = len(j)
                    max_twoc = len(j)
        for j in i:
            if len(j) < max_twoc:
                if len(j) > max_threec:
                    crystalthree = len(j)
                    max_threec = len(j)
        print 'frame is', k
        print 'max crystal size is crystal', len(crystal)
        print 'crystal two',crystaltwo
        crystal_cords_frame.append(crystal_cords)
        max_crystal.append(len(crystal))
        second_crystal.append(crystaltwo)
        third_crystal.append(crystalthree)
        num_crystal.append(count)


    f = open('large.xyz','w')
    count=0
    print 'crystal cord frame'
    print len(crystal_cords_frame)
    print M.shape
    count = 0
    for k in range(0,M.shape[0]):
        print k
        f.write(('%i \n\n')%(M.shape[1]))
        for i in range(M.shape[1]):
            if i in crystal_cords_frame[k]:
                if M[k][i][3] == 0:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',M[k][i][0],M[k][i][1],M[k][i][2]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',M[k][i][0],M[k][i][1],M[k][i][2]))
            else:
                L = 80
                if M[k][i][3] == 0:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',4*L[0],4*L[0],4*L[0]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',4*L[0],4*L[0],4*L[0]))
    f.close()



    x = range(num_frames)
    xlabel = 'frames'
    ylabel = 'sigma'
    util.pickle_dump([num_crystal,max_crystal,second_crystal,third_crystal],'numsizes.pkl')
    pyplot.plot3(x,num_crystal,x,max_crystal,x,second_crystal,xlabel=xlabel,ylabel=ylabel,label1='num',label2='max',label3='second',showleg=True,save='crystal_sizes')
#Find the number of solid particles within the system
def solid_particles(VW,V,W,L,VW_names,skip=30,step=5e4,rcut=False):
    BCC = []
    SC = []
    SCW = []
    BCCcrystals = []
    SCcrystals = []
    SCWcrystals = []
    try:
        SCcrystals = util.pickle_load('sccrystal.pkl')
        SCWcrystals = util.pickle_load('scWcrystal.pkl')
        BCCcrystals = util.pickle_load('bccrystal.pkl')
        SC = util.pickle_load('sc.pkl')
        SCW = util.pickle_load('scW.pkl')
        BCC = util.pickle_load('bcc.pkl')
    except:
        if rcut == False:
            rcutVW=min_particle_distance(VW[-1:],VW[-1:],L)+3
        rcutV = rcutVW*2/3**0.5
        print 'rcut'
        #sp 8
        #flexible
        rcutVW = 17
        rcutV = 20
        print rcutV
        print rcutVW

        #rigid
        #rcutVW = 23
        #rcutV = 28
        for i in range(0,VW.shape[0],skip):
            print 'step', i
            avg = 5
            num, crystals = (bond_order.solid_particles(av.average_position(VW[i:i+avg]),
                            L,c_cut=0.1, rcut=rcutVW,l=6,crystal=6))
            BCCcrystals.append(crystals)
            BCC.append(num)
            num, crystals = (bond_order.solid_particles(av.average_position(V[i:i+avg]),
                L,rcut=rcutV,l=4,  c_cut=0.1,crystal=3))
            numW, crystalsW =  (bond_order.solid_particles(av.average_position(W[i:i+avg]),
                L,rcut=rcutV,l=4,  c_cut=0.1,crystal=3))
            SC.append(num)
            SCcrystals.append(crystals)
            SCW.append(numW)
            SCWcrystals.append(crystalsW)
        util.pickle_dump(SCcrystals,'sccrystal.pkl')
        util.pickle_dump(SCWcrystals,'scWcrystal.pkl')
        util.pickle_dump(BCCcrystals,'bccrystal.pkl')
        util.pickle_dump(SC,'sc.pkl')
        util.pickle_dump(SCW,'scW.pkl')
        util.pickle_dump(BCC,'bcc.pkl')

    x = np.arange(len(BCC))
    x *= skip
    xlabel ='time'
    ylabel ='solid'
    SC = np.array(SC)/float(V.shape[1])
    SCW = np.array(SCW)/float(W.shape[1])
    BCC = np.array(BCC)/float(VW.shape[1])
    pyplot.plot3(x,SC,x,SCW,x,BCC, xlabel=xlabel,ylabel=ylabel,save='solidcrystal')
    # make an xyz file of the crystals and there locations
    f = open('bcccrystals.xyz','w')
    count=0
    for k in range(0,VW.shape[0],skip):
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(BCCcrystals[count]))):
            if BCCcrystals[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()
    f = open('sccrystals.xyz','w')
    count=0
    for i in range(0,VW.shape[0],skip):
        f.write(('%i \n')%(len(SCcrystals[count])))
        f.write('Atoms\n')
        for j in range(len(SCcrystals[count])):
            if SCcrystals[count][j] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()

    f = open('scWcrystals.xyz','w')
    count=0
    for i in range(0,VW.shape[0],skip):
        f.write(('%i \n')%(len(SCWcrystals[count])))
        f.write('Atoms\n')
        for j in range(len(SCWcrystals[count])):
            if SCWcrystals[count][j] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()
    return x,BCC,SC
#Find the number of solid particles within the system
def solid_particles_VW(VW,V,W,L,VW_names,skip=30,step=5e4,rcut=False):
    BCC = []
    BCCcrystals = []
    try:
        BCCcrystals = util.pickle_load('bccrystal.pkl')
        BCC = util.pickle_load('bcc.pkl')
    except:
        if rcut == False:
            rcutVW=min_particle_distance(VW[-1:],VW[-1:],L)+3
        rcutV = rcutVW*2/3**0.5
        print 'rcut'
        #sp 8
        #flexible
        rcutVW = 20
        print rcutV
        print rcutVW

        #rigid
        #rcutVW = 23
        #rcutV = 28
        for i in range(0,VW.shape[0],skip):
            print 'step', i
            avg = 5
            #num, crystals = (bond_order.solid_particles(av.average_position(VW[i:i+avg]),
            #                L,c_cut=0.1, rcut=rcutVW,l=6,crystal=6))
            num, crystals = (bond_order.solid_particles(VW[i],
                            L,c_cut=0.1, rcut=rcutVW,l=6,crystal=6))
            BCCcrystals.append(crystals)
            BCC.append(num)
        util.pickle_dump(BCCcrystals,'bccrystal.pkl')
        util.pickle_dump(BCC,'bcc.pkl')

    x = np.arange(len(BCC))
    x *= skip
    xlabel ='time'
    ylabel ='solid'
    BCC = np.array(BCC)/float(VW.shape[1])
    pyplot.plot(x,BCC, xlabel=xlabel,ylabel=ylabel,save='solidcrystal')
    # make an xyz file of the crystals and there locations
    f = open('bcccrystals.xyz','w')
    count=0
    for k in range(0,VW.shape[0],skip):
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(BCCcrystals[count]))):
            if BCCcrystals[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()
    return x,BCC
def bcc_particles_VW(VW,V,W,L,VW_names,skip=30,step=5e4,rcut=False):
    BCC = []
    BCCcrystals = []
    try:
        BCCcrystals = util.pickle_load('bccrystal.pkl')
        BCC = util.pickle_load('bcc.pkl')
    except:
        if rcut == False:
            rcutVW=min_particle_distance(VW[-1:],VW[-1:],L)+3
        rcutV = rcutVW*2/3**0.5
        print 'rcut'
        #sp 8
        #flexible
        rcutVW = 20
        print rcutV
        print rcutVW

        #rigid
        #rcutVW = 23
        #rcutV = 28
        for i in range(0,VW.shape[0],skip):
            print 'step', i
            avg = 5
            #num, crystals = (bond_order.solid_particles(av.average_position(VW[i:i+avg]),
            #                L,c_cut=0.1, rcut=rcutVW,l=6,crystal=6))
            num, crystals = (bond_order.solid_particles(VW[i],
                            L,c_cut=0.1, rcut=rcutVW,l=6,crystal=6))
            BCCcrystals.append(crystals)
            BCC.append(num)
        util.pickle_dump(BCCcrystals,'bccrystal.pkl')
        util.pickle_dump(BCC,'bcc.pkl')

    x = np.arange(len(BCC))
    x *= skip
    xlabel ='time'
    ylabel ='solid'
    BCC = np.array(BCC)/float(VW.shape[1])
    pyplot.plot(x,BCC, xlabel=xlabel,ylabel=ylabel,save='solidcrystal')
    # make an xyz file of the crystals and there locations
    f = open('bcccrystals.xyz','w')
    count=0
    for k in range(0,VW.shape[0],skip):
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(BCCcrystals[count]))):
            if BCCcrystals[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()
    return x,BCC
#Find the number of solid particles within the system
def solid_particles_surround(VW,V,W,L,VW_names,skip=30,step=5e4,rcut=False):
    try:
        surface = util.pickle_load('surface.pkl')
    except:
        ###################################
        # normal crystal
        #BCCcrystals = util.pickle_load('bccrystal.pkl')
        #BCC = util.pickle_load('bcc.pkl')
        #####################################
        ### filtered crystal
        from MD.util import readxyz
        fid = open('large.xyz','r')
        M = readxyz.ReadCord(trajectory = 'large.xyz',frames = 777)
        crystal = M.cord_auto(['V','W'])
        bcc = np.zeros((777,432,1))
        for frame in range(bcc.shape[0]):
            for i in range(bcc.shape[1]):
                if crystal[frame][i][0] > L[0]:
                    bcc[frame][i]=0
                else:
                    bcc[frame][i]=1
        ######################################
        BCCcrystals = bcc
        surface = np.zeros((len(BCCcrystals),len(BCCcrystals[0])))
        for k in range(len(BCCcrystals)):
            print k
            for i in range(len(BCCcrystals[k])):
                #if 4 neighbors are crystals, we call it the surface particle 
                if BCCcrystals[k][i] == 0:
                    nn = nearest_neighbors_index(VW[k],i,L,rcut=17)[0]
                    count = 0
                    for j in nn:
                        if BCCcrystals[k][j] == 1:
                            count += 1
                    if count >=2:
                        surface[k][i] = 1
        util.pickle_dump(surface,'surface.pkl')

    ###############
    # write out surface
    ################
    x = np.arange(surface.shape[0])
    x *= skip
    xlabel ='time'
    ylabel ='solid'
    # make an xyz file of the crystals and there locations
    f = open('surfacecrystals.xyz','w')
    count=0
    print len(surface)
    for k in range(0,len(surface),skip):
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(surface[count]))):
            if surface[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()
    return surface
###########################################################
# solid particles specifiallcy for cubes
#########################################################
def solid_particles_sc(VW,V,W,L,skip=30,step=5e4,rcut=False):
    SC = []
    SCcrystals = []
    try:
        SCcrystals = util.pickle_load('bccrystal.pkl')
        SC = util.pickle_load('sc.pkl')
    except:
        #cutoff for a crystal
        if rcut == False:
            rcutVW=min_particle_distance(VW[-1:],VW[-1:],L)+3
        rcutV = rcutVW*2/3**0.5
        rcutVW = 17
        rcutV = 20

        for i in range(0,VW.shape[0],skip):
            print 'step', i
            num, crystals = (bond_order.solid_particles(VW[i],L,c_cut=0.1,
                rcut=rcutVW,l=4,crystal=3))
            SCcrystals.append(crystals)
            SC.append(num)
        util.pickle_dump(SCcrystals,'sccrystal.pkl')
        util.pickle_dump(SC,'sc.pkl')

    x = np.arange(len(BCC))
    x *= skip
    xlabel ='time'
    ylabel ='solid'
    SC = np.array(SC)/float(V.shape[1])
    SCW = np.array(SCW)/float(W.shape[1])
    BCC = np.array(BCC)/float(VW.shape[1])
    pyplot.plot3(x,SC,x,SCW,x,BCC, xlabel=xlabel,ylabel=ylabel,save='solidcrystal')
    # make an xyz file of the crystals and there locations
    f = open('sccrystals.xyz','w')
    count=0
    for i in range(0,VW.shape[0],skip):
        f.write(('%i \n')%(len(SCcrystals[count])))
        f.write('Atoms\n')
        for j in range(len(SCcrystals[count])):
            if SCcrystals[count][j] == 1:
                if j<VW.shape[1]/2:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',VW[i][j][0],VW[i][j][1],VW[i][j][2]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',VW[i][j][0],VW[i][j][1],VW[i][j][2]))
            else:
                if j<VW.shape[1]/2:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',4*L[0],4*L[0],4*L[0]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',4*L[0],4*L[0],4*L[0]))
        count += 1
    f.close()

    return x,SC
#Find the number of crystals, the size of the three largest crystals, and the width of the
#largest crystal in the x, y, and z direction
def crystal_size_cube(L,num_frames,num_particles,rcut=15):
    #load all crystal positions
    #num_frames = 777
    #rigid
    rcut = 25
    fid = open('sccrystals.xyz','r')
    M = np.zeros((num_frames,num_particles,4))
    for frame in range(num_frames):
        fid.readline()
        fid.readline()
        for i in range(num_particles):
            row = fid.readline().split()
            M[frame][i][0]=float(row[1])
            M[frame][i][1]=float(row[2])
            M[frame][i][2]=float(row[3])
            if row[0] == 'V':
                M[frame][i][3]=0
            if row[0] == 'W':
                M[frame][i][3]=1
    #Find conected crystals 
    import MD.analysis.connections as con
    import MD.analysis.graph as graph
    import networkx as nx
    try:
        networks = util.pickle_load('cnetwork.pkl')
    except:
        connections = con.connections_min_max(M,M,L,rcut=rcut,rmin=4)
        for i in range(len(connections)):
            for j in range(len(connections[i])):
                connections[i][j][1] -= num_particles
        #Find groups of crystals via percolation, ie use distance as connection
        networks, num_networks, deg, neighbors, num_n, gr =   graph.grapher(connections,num_particles,1)
        util.pickle_dump(networks,'cnetwork.pkl')
    #Once groups of crystals are identified, find the size of each crystal
    #TO-DO
    #find the largest crystal in each frame
    num_crystal = []
    max_crystal = []
    second_crystal = []
    third_crystal = []
    crystal_cords_frame = []
    for k,i in enumerate(networks):
        crystal = []
        crystaltwo = 0
        crystalthree = 0
        max_c = 1
        max_twoc = 1
        max_threec = 1
        count = 0
        for j in i:
            if len(j) > max_c:
                crystal = j
                max_c = len(j)
                crystal_cords = j
            if len(j) > 1:
                count += 1
        for j in i:
            if len(j) < max_c:
                if len(j) > max_twoc:
                    crystaltwo = len(j)
                    max_twoc = len(j)
        for j in i:
            if len(j) < max_twoc:
                if len(j) > max_threec:
                    crystalthree = len(j)
                    max_threec = len(j)
        print 'frame is', k
        print 'max crystal size is crystal', len(crystal)
        print 'crystal two',crystaltwo
        crystal_cords_frame.append(crystal_cords)
        max_crystal.append(len(crystal))
        second_crystal.append(crystaltwo)
        third_crystal.append(crystalthree)
        num_crystal.append(count)


    f = open('large.xyz','w')
    count=0
    print 'crystal cord frame'
    print len(crystal_cords_frame)
    print M.shape
    count = 0
    for k in range(0,M.shape[0]):
        print k
        f.write(('%i \n\n')%(M.shape[1]))
        for i in range(M.shape[1]):
            if i in crystal_cords_frame[k]:
                if M[k][i][3] == 0:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',M[k][i][0],M[k][i][1],M[k][i][2]))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',M[k][i][0],M[k][i][1],M[k][i][2]))
            else:
                L = 80
                if M[k][i][3] == 0:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',4*L,4*L,4*L))
                else:
                    f.write(('%c %.2f %.2f %.2f \n')%('W',4*L,4*L,4*L))
    f.close()



    x = range(num_frames)
    xlabel = 'frames'
    ylabel = 'sigma'
    util.pickle_dump([num_crystal,max_crystal,second_crystal,third_crystal],'numsizes.pkl')
    pyplot.plot3(x,num_crystal,x,max_crystal,x,second_crystal,xlabel=xlabel,ylabel=ylabel,label1='num',label2='max',label3='second',showleg=True,save='crystal_sizes')



if __name__ == '__main__':
    #run_debug()
    #run_all()
    import MD
    reload(bond_order)
    solid_particles(MD.VW,MD.V,MD.W,MD.L,skip=2,rcut=False)
    crystal_size(L,VW.shape[0],VW.shape[1])
