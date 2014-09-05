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
# solid particles specifiallcy for bcc systems
#########################################################
#Find the number of crystals, the size of the three largest crystals, and the width of the
#largest crystal in the x, y, and z direction
def Average(VW,n_start,n_frames,L):
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
    fid = open('average.xyz','w')
    fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
    for i in range(Average.shape[0]):
        fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
    fid.close()
    return Average
    
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

def q6q4_map(VW,L,x,avg=1,step=5e4,rcut=30):
    Q4 = []
    Q6 = []
    xlabel ='Q4'
    ylabel ='Q6'
    try:
        Q6 = util.pickle_load('Q6.pkl')
        Q4 = util.pickle_load('Q4.pkl')
        print 'plotting'
        for i in range(len(x)):
            pyplot.plot(Q4[i],Q6[i],
                    xlabel=xlabel,ylabel=ylabel,save='Qmap_frame%i'%x[i],
                    limx=(0,.2),limy=(0,.6),linestyle='',marker='x')
    except:
        for k in x:
            print 'step', k
            if avg > 1:
                A = Average(VW,k,avg,L)
            else:
                A = VW[k]
            Q4_avg = (bond_order.Qmap(A,L[k],count=14,l=4,rcut=rcut))
            Q6_avg = (bond_order.Qmap(A,L[k],count=14,l=6,rcut=rcut))
            Q4.append(Q4_avg.real)
            Q6.append(Q6_avg.real)
            pyplot.plot(Q4_avg,Q6_avg,
                    xlabel=xlabel,ylabel=ylabel,save='Qmap_frame%i'%k,
                    limx=(0,.2),limy=(0,.6),linestyle='',marker='x')
    fid = open('Q6Q4map.txt','w')
    for i in range(len(x)):
        fid.write('Step %i length %.2f'%(x[i],L[i][0]))
        for j in range(Q6[i].shape[0]):
            fid.write('%i %.8f %.8f\n'%(j,Q6[i][j],Q4[i][j]))
    fid.close()
    return Q4,Q6
    
def write_qmap_xyz(VW,L,x,Q6,Q4):
    fid = open('q6q4.xyz','w')
    bcc_cut = .05
    bcc = np.array([.45,.03 ])
    hcp_cut = .035
    hcp = np.array([.475,.1 ])
    fcc_cut = .075
    fcc = np.array([.5,.17 ])
    Liq_cut = .1
    Liq = np.array([.2,.04])
    C = np.array([200,200,200])
    for i in range(len(x)):
        fid.write('%i\n length %.2f\n'%(5*VW.shape[1],L[x[i]][0]))
        Qx = [[]  for n in range(5)]
        Qy = [[]  for n in range(5)]
        for j in range(VW.shape[1]):
            w = [False for n in range(5)]
            #bcc
            if points.dist_np(bcc,np.array([Q6[i][j],Q4[i][j]]))[0] < bcc_cut:
                w[0] = True
            #hcp
            elif points.dist_np(hcp,np.array([Q6[i][j],Q4[i][j]]))[0] < hcp_cut:
                w[1] = True
            #fcc
            elif points.dist_np(fcc,np.array([Q6[i][j],Q4[i][j]]))[0] < fcc_cut:
                w[2] = True
            #liq
            elif points.dist_np(Liq,np.array([Q6[i][j],Q4[i][j]]))[0] < Liq_cut:
                w[3] = True
            #Interface
            else:
                w[4] = True
            if w[0] == True:
                fid.write(('B %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
                    VW[x[i]][j][1],VW[x[i]][j][2])))
                Qx[0].append(Q4[i][j])
                Qy[0].append(Q6[i][j])
            else:
                fid.write(('B %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
            if w[1] == True:
                fid.write(('H %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
                    VW[x[i]][j][1],VW[x[i]][j][2])))
                Qx[1].append(Q4[i][j])
                Qy[1].append(Q6[i][j])
            else:
                fid.write(('H %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
            if w[2] == True:
                fid.write(('F %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
                    VW[x[i]][j][1],VW[x[i]][j][2])))
                Qx[2].append(Q4[i][j])
                Qy[2].append(Q6[i][j])
            else:
                fid.write(('F %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
            if w[3] == True:
                fid.write(('L %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
                    VW[x[i]][j][1],VW[x[i]][j][2])))
                Qx[3].append(Q4[i][j])
                Qy[3].append(Q6[i][j])
            else:
                fid.write(('L %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
            if w[4] == True:
                fid.write(('I %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
                    VW[x[i]][j][1],VW[x[i]][j][2])))
                Qx[4].append(Q4[i][j])
                Qy[4].append(Q6[i][j])
            else:
                fid.write(('I %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
        print 'plotting'
        xlabel = 'q4_avg'
        ylabel = 'q6_avg'
        Label = ['bcc','hcp','fcc','liq','int']
        print x[i]
        pyplot.plot_multi(Qx,Qy,Label,
                xlabel=xlabel,ylabel=ylabel,save='Qmap_frame%i'%x[i],
                limx=(0,.2),limy=(0,.6),make_marker=True)
    fid.close()
def q6q4_fast_map(VW,L,x,avg=1,step=5e4,rcut=30):
    Q4 = []
    Q6 = []
    xlabel ='Q4'
    ylabel ='Q6'
    try:
        Q6 = util.pickle_load('Q6.pkl')
        Q4 = util.pickle_load('Q4.pkl')
        #print 'plotting'
        #for i in range(len(x)):
        #    pyplot.plot(Q4[i],Q6[i],
        #            xlabel=xlabel,ylabel=ylabel,save='Qmap_frame%i'%x[i],
        #            limx=(0,.2),limy=(0,.6),linestyle='',marker='x')
    except:
        for k in x:
            print 'step', k
            if avg > 1:
                A = Average(VW,k,avg,L)
            else:
                A = VW[k]
            Q4_avg,Q6_avg = (bond_order.QmapFast(A,L[k],count=24,rcut=rcut))
            Q4.append(Q4_avg.real)
            Q6.append(Q6_avg.real)
            pyplot.plot(Q4_avg,Q6_avg,
                    xlabel=xlabel,ylabel=ylabel,save='Qmap_frame%i'%k,
                    limx=(0,.2),limy=(0,.6),linestyle='',marker='x')
        util.pickle_dump(Q6,'Q6.pkl')
        util.pickle_dump(Q4,'Q4.pkl')
    fid = open('Q6Q4map.txt','w')
    for i in range(len(x)):
        fid.write('Step %i length %.2f'%(x[i],L[i][0]))
        for j in range(Q6[i].shape[0]):
            fid.write('%i %.8f %.8f\n'%(j,Q6[i][j],Q4[i][j]))
    fid.close()
    #fid = open('q6q4.xyz','w')
    #bcc_cut = .05
    #bcc = np.array([.45,.03 ])
    #hcp_cut = .035
    #hcp = np.array([.475,.1 ])
    #fcc_cut = .075
    #fcc = np.array([.5,.17 ])
    #Liq_cut = .1
    #Liq = np.array([.2,.04])
    #for i in range(len(x)):
    #    fid.write('%i\n length %.2f\n'%(VW.shape[1],L[x[i]][0]))
    #    for j in range(VW.shape[1]):
    #        #bcc
    #        if points.dist_np(bcc,np.array([Q6[i][j],Q4[i][j]]))[0] < bcc_cut:
    #            fid.write(('B %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
    #                VW[x[i]][j][1],VW[x[i]][j][2])))
    #        #hcp
    #        elif points.dist_np(hcp,np.array([Q6[i][j],Q4[i][j]]))[0] < hcp_cut:
    #            fid.write(('H %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
    #                VW[x[i]][j][1],VW[x[i]][j][2])))
    #        #fcc
    #        elif points.dist_np(fcc,np.array([Q6[i][j],Q4[i][j]]))[0] < fcc_cut:
    #            fid.write(('F %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
    #                VW[x[i]][j][1],VW[x[i]][j][2])))
    #        #liq
    #        elif points.dist_np(Liq,np.array([Q6[i][j],Q4[i][j]]))[0] < Liq_cut:
    #            fid.write(('L %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
    #                VW[x[i]][j][1],VW[x[i]][j][2])))
    #        #Interface
    #        else:
    #            fid.write(('I %.4f %.4f %.4f\n'%(VW[x[i]][j][0],
    #                VW[x[i]][j][1],VW[x[i]][j][2])))
    #fid.close()
    print 'writing'
    write_qmap_xyz(VW,L,x,Q6,Q4)
    print 'finished writing'
    return Q4,Q6
def steinhardt(VW,L,x,avg=1,step=5e4,rcut=30):
    Q4 = []
    Q6 = []
    for k in x:
        print 'step', k
        if avg > 1:
            A = Average(VW,k,avg,L)
        else:
            A = VW[k]
        Q2_avg = (bond_order.steinhardt_order(A,L[k],l=2,rcut=rcut))
        Q4_avg = (bond_order.steinhardt_order(A,L[k],l=4,rcut=rcut))
        Q6_avg = (bond_order.steinhardt_order(A,L[k],l=6,rcut=rcut))
        Q8_avg = (bond_order.steinhardt_order(A,L[k],l=8,rcut=rcut))
        Q10_avg = (bond_order.steinhardt_order(A,L[k],l=10,rcut=rcut))
#Find the number of solid particles within the system
def solid_particles_VW(VW,L,VW_names,x,avg=1,step=5e4):
    BCC = []
    BCCcrystals = []
    try:
        BCCcrystals = util.pickle_load('bccrystal.pkl')
        BCC = util.pickle_load('bcc.pkl')
        asdasdf
    except:
        for k in x:
            print 'step', k
            A = Average(VW,k,avg,L)
            num, crystals = (bond_order.solid_particles(A,
                            L[k],c_cut=0.15, count=8,l=6,crystal=7))
            BCCcrystals.append(crystals)
            BCC.append(num)
        util.pickle_dump(BCCcrystals,'bccrystal.pkl')
        util.pickle_dump(BCC,'bcc.pkl')

    fid = open('solidpart.txt','w')
    for i in range(len(x)):
        fid.write('%i %i %.2f\n'%(x[i],BCC[i],L[x[i]][0]))
    fid.close()
    xlabel ='time'
    ylabel ='solid'
    BCC = np.array(BCC)/float(VW.shape[1])
    pyplot.plot(x,BCC, xlabel=xlabel,ylabel=ylabel,save='solidcrystal')
    # make an xyz file of the crystals and there locations
    f = open('bcccrystals.xyz','w')
    count=0
    for k in x:
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(BCCcrystals[count]))):
            if BCCcrystals[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f\n')%('W',4*L[k][0],4*L[k][0],4*L[k][0]))
        count += 1
    f.close()
    return x,BCC
#Find the number of solid particles within the system
#Find the number of solid particles within the system
def solid_particles_sc(VW,L,VW_names,x,step=5e4):
    SC = []
    SCcrystals = []
    try:
        SCcrystals = util.pickle_load('sccrystal.pkl')
        SC = util.pickle_load('sc.pkl')
    except:
        for i in range(len(x)-1):
            a = x[i+1]-x[i]
            delta =  a/4
            k = x[i]+a/2
            print delta
        for k in x:
            print 'step', k
            num, crystals = (bond_order.solid_particles(VW[k],
                            L[k],c_cut=0.175, count=6,l=4,crystal=5))
            SCcrystals.append(crystals)
            SC.append(num)
        util.pickle_dump(SCcrystals,'sccrystal.pkl')
        util.pickle_dump(SC,'sc.pkl')

    fid = open('solidpartsc.txt','w')
    for i in range(len(x)):
        fid.write('%i %i\n'%(x[i],SC[i]))
    fid.close()
    xlabel ='time'
    ylabel ='solid'
    SC = np.array(SC)/float(VW.shape[1])
    pyplot.plot(x,SC, xlabel=xlabel,ylabel=ylabel,save='solidcrystalsc')
    # make an xyz file of the crystals and there locations
    f = open('sccrystals.xyz','w')
    count=0
    for k in x:
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(SCcrystals[count]))):
            if SCcrystals[count][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%('V',VW[k][i][0],VW[k][i][1],VW[k][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f\n')%('W',4*L[k][0],4*L[k][0],4*L[k][0]))
        count += 1
    f.close()
    return x,SC
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
def sc_surface(VW,L,rcut=15,delta=25,c_cut=2):
    #####################################
    c_index = util.pickle_load('sccrystal.pkl')
    surface = np.zeros((len(c_index),len(c_index[0])))
    ######################################
    for k in range(len(c_index)):
        for i in range(len(c_index[k])):
            #if c_cut neighbors are crystals, we call it the surface particle 
            if c_index[k][i] == 0:
                print k
                print k*delta
                print VW.shape

                nn = nearest_neighbors_index(VW[k*delta],i,L[k*delta],rcut=rcut)[0]
                count = 0
                for j in nn:
                    if c_index[k][j] == 1:
                        count += 1
                if count >=c_cut:
                    surface[k][i] = 1
    util.pickle_dump(surface,'surface.pkl')
    counter = []
    fid = open('surfacecount.txt','w')
    for k in range(surface.shape[0]):
        count = 0
        for i in range(surface.shape[1]):
            if surface[k][i] == 1:
                count+=1
        fid.write('%i %i\n'%(k,count))

    ###############
    # write out surface
    ################
    xlabel ='time'
    ylabel ='solid'
    # make an xyz file of the crystals and there locations
    f = open('surfacecrystals.xyz','w')
    print len(surface)
    VW_names  = []
    for i in range(VW.shape[1]):
        if i <VW.shape[1]/2.:
            VW_names.append('V')
        else:
            VW_names.append('W')
    for k in range(0,len(surface)):
        f.write(('%i \n\n')%(VW.shape[1]))
        for i in range((len(surface[k]))):
            if surface[k][i] == 1:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],
                        VW[k*delta][i][0],VW[k*delta][i][1],VW[k*delta][i][2]))
            else:
                    f.write(('%c %.2f %.2f %.2f \n')%(VW_names[i],
                        4*L[0][0],4*L[0][0],4*L[0][0]))
    f.close()
    return surface
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
    import MD.unit.make_bcc as make_bcc
    a = 17.
    basis = np.array([0.0,0.0,0.0])
    a1 = np.array([-0.5*a,0.5*a,0.5*a])
    a2 = np.array([0.5*a,-0.5*a,0.5*a])
    a3 = np.array([0.5*a,0.5*a,-0.5*a])
    L = [np.array([a*5,a*5,a*5])]
    VW, VW_names = make_bcc.make_bcc(a1,a2,a3,basis,L[0],S=30,name='bcc_unit_large.xyz')
    reload(bond_order)
    q6q4_fast_map(np.array([VW]),L,[0],rcut=30)
    steinhardt(np.array([VW]),L,[0],rcut=30)
    #crystal_size(L,VW.shape[0],VW.shape[1])
