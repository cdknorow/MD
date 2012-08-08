import numpy as np
import MD.util as util
from util import readxyz
import MD.base.points as points
import MD.plot.pyplot as pyplot
import matplotlib.pyplot as plt
import MD.unit.make_bcc as make
import MD.analysis.sfactor as sf
import MD.analysis.msdr as msdr
import MD.plot.pyplot as pyplot
import MD.analysis.nearest_neighbor as nn
from MD.analysis.particle_distance import min_point_distance
reload(nn)
from MD.analysis.nearest_neighbor import nearest_neighbors_point
from MD.analysis.nearest_neighbor import close_neighbors_point
from MD.analysis.nearest_neighbor import rcut_neighbors_point
from MD.plot.histogram import histogram_reg 
# find the defects in the lattice
#
# Returns DV, DW, the types of defects at each frame 1, 0, 2 (none, vacancy,
# interstitial) - represent substitutional.
#
#for non comensurate systems
def find_defects_mod(CV,CW,VW,V,W,L,n_finish=1,n_start=0,delta=20,filters=0.2):
    from MD.analysis.nearest_neighbor import ws_neighbors_point
    from MD.analysis.nearest_neighbor import close_neighbors_point
    #######
    debug = False
    x = np.arange(n_start,n_finish,delta)
    lattice = []
    for i in range(CV.shape[0]):
        lattice.append(int(points.dist(CV[0],CV[i],L)[0]))
    for i in range(CW.shape[0]):
        lattice.append(int(points.dist(CV[0],CW[i],L)[0]))
    points.unique(lattice)
    lattice.sort()
    print lattice
    rmin = lattice[1] / 2.0
    rmax = lattice[2] / 2.0 +1
    #for interstitial
    s = min(W.shape[1],V.shape[1])
    #for vacancy
    #s = max(W.shape[1],V.shape[1])
    DW = np.zeros((len(x),CW.shape[0]))
    DV = np.zeros((len(x),CV.shape[0]))
    vac_count = np.zeros((len(x),1))
    int_count = np.zeros((len(x),1))
    fid = open('defects.txt','w')
    master = [[],[]]
    particle_frame = []
    for index,k in enumerate(x):
        print 'frame',k
        # Identify the points on the bcc Lattice which belong to A and B type
        #lets look at one point first
        # Assign each point to its place on the lattice
        #find the closest point   
        #####################
        N = []
        N_V = []
        N_W = []
        Vn = []
        Wn = []
        for i in V[k]:
            if abs(i[0]) < L[0]:
                Vn.append(i)
        for i in W[k]:
            if abs(i[0]) < L[0]:
                Wn.append(i)
        print 'Wn',len(Wn)
        print 'Vn',len(Vn)
        particle_frame.append(len(Wn)+len(Vn))
        Wn = np.array(Wn)
        Vn = np.array(Vn)

        for i in range(CV.shape[0]):
            N =  ws_neighbors_point(Vn,CV[i],L,i,rmin=rmin,rmax=rmax)[0]
            N_V.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(Wn,CV[i],L,i,rmin=rmin,rmax=rmax)[0]
                N_W.extend(N2)
                num = -len(N2)
            DV[index][i] = num
        ###########################
        for i in range(CW.shape[0]):
            N =  ws_neighbors_point(Wn,CW[i],L,i+W.shape[1],rmin=rmin,rmax=rmax)[0]
            N_W.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(Vn,CW[i],L,i+V.shape[1],rmin=rmin,rmax=rmax)[0]
                N_V.extend(N2)
                num = -len(N2)
            DW[index][i] = num
        #find the atoms that haven't been placed on the lattice yet
        IV = points.difference(range(Vn.shape[0]),N_V)
        IW = points.difference(range(Wn.shape[0]),N_W)
        ##
        if debug:
            print 'atoms not added listed after first sorting'
            print IV, IW
            print DW[index]
            print DV[index]
        for i in IV:
            #find closest lattice point
            nw, dw = close_neighbors_point(CW,Vn[i],L)
            nv, dv = close_neighbors_point(CV,Vn[i],L)
            if dw <= dv:
                #check to see if there is already an atom at that point
                if DW[index][nw] == 1 or DW[index][nw] == -1:
                    if DV[index][nv] == 0:
                        DV[index][nv] = 1
                        N_V.extend([i])
                    else:
                        DW[index][nw] = -2
                        N_V.extend([i])
                if DW[index][nw] == 0:
                        DW[index][nw] = -1
                        N_V.extend([i])
                #check to see if there is already an atom at that point
            else:
                if DV[index][nv] == 1 or DV[index][nv] == -1:
                    #if there isn't one at the other point add it
                    if DW[index][nw] == 0:
                        DW[index][nw] = -1
                        N_V.extend([i])
                    else:
                        if DV[index][nv] == 1:
                            DV[index][nv] = 2
                        if DV[index][nv] == -1:
                            DV[index][nv] = -2
                        N_V.extend([i])
                if DV[index][nv] == 0:
                    DV[index][nv] = 1
                    N_V.extend([i])
        for i in IW:
            nw, dw = close_neighbors_point(CW,Wn[i],L)
            nv, dv = close_neighbors_point(CV,Wn[i],L)
            if dv <= dw:
                if DV[index][nv] == 1 or DV[index][nv] == -1:
                    if DW[index][nw] == 0:
                        DW[index][nw] = 1
                        N_W.extend([i])
                    else:
                        DV[index][nv] = -2
                        N_W.extend([i])
                if DV[index][nv] == 0:
                    DV[index][nv] = -1
                    N_W.extend([i])
            else:
                if DW[index][nw] == 1 or DW[index][nw] == -1:
                    if DV[index][nv] == 0:
                        DV[index][nv] = -1
                        N_W.extend([i])
                    else:
                        DW[index][nw] = 2
                        N_W.extend([i])
                if DW[index][nw] == 0:
                    DW[index][nw] = 1
                    N_W.extend([i])
        #find the atoms that haven't been placed on the lattice yet
        IV = points.difference(range(Vn.shape[0]),N_V)
        IW = points.difference(range(Wn.shape[0]),N_W)
        ##
        if debug:
            print 'atoms not added list for debugging'
            print IV, IW
            print DW[index]
            print DV[index]
            print 'Defect list for lattice'
        #print out the vacency, substitutions
        def out_defect(A, index, fid, def_list, C=0):
            for i in range(A.shape[1]):
                if A[index][i] == 0:
                    pr =  'vacecy at '+ str(i+C)+ '\n'
                    try:
                        def_list[0].extend(i+C)
                    except:
                        def_list[0].append(i+C)
                    if debug:
                        print pr
                    fid.write(pr)
                if A[index][i] == -1:
                    pr = 'substitution ' + str(i+C)+ '\n'
                    if debug:
                        print pr
                    fid.write(pr)
                if A[index][i] == -2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
                    if debug:
                        print pr
                    fid.write(pr)
                if A[index][i] == 2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
                    if debug:
                        print pr
                    fid.write(pr)

        frame = 'Frame ' + str(k) + '\n'
        fid.write(frame)
        def_list = [[],[]]
        out_defect(DV, index, fid, def_list)
        out_defect(DW, index, fid, def_list, C = DV.shape[1])
        if len(points.difference(def_list[0], master[0])) != 0:
            vac_count[index] += len(points.difference(def_list[1],master[1]))
            master[0] = def_list[0]
        if len(points.difference(def_list[1],master[1])) != 0:
            int_count[index] += len(points.difference(def_list[1],master[1]))
            master[1] = def_list[1]

        #find the atoms that haven't been placed on the lattice yet
        IV = points.difference(range(V.shape[1]),N_V)
        IW = points.difference(range(W.shape[1]),N_W)

    # Identify the  defects surrounding each point
    #The total number of frames we are going to look at
    # Find the number of defects in each frame
    count = 0
    substitutions = []
    vacancies = []
    intersticial = []
    def count_def(A,check):
        count = 0
        for i in A:
            if i == check:
                count +=1
        return count
    def count_ldef(A,check):
        count = 0
        for i in A:
            if i < check:
                count +=1
        return count
    def count_adef(A,check):
        count = 0
        for i in A:
            if abs(i) == check:
                count +=1
        return count
    for k in range(DW.shape[0]):
        #find the points that are nearest neighbor that are different
        substitutions.append(count_ldef(DW[k],0)+count_ldef(DV[k],0))
        vacancies.append(count_def(DW[k],0)+count_def(DV[k],0))
        intersticial.append(count_adef(DW[k],2.0)+count_adef(DV[k],2.0))
    print 'substitions'
    print substitutions
    print 'vacancies'
    print vacancies
    print 'interstitials'
    print intersticial
    util.pickle_dump(DV,'DV_s.pkl')
    util.pickle_dump(DW,'DW_s.pkl')

    util.pickle_dump([substitutions, vacancies, intersticial, x],'plot_def_s.pkl')
    pyplot.plot3(x, substitutions, x, particle_frame, x, intersticial,
            label1='substitutions',label2='number of particles',label3='intersticial',
            save='defects_per_frame_surface',showleg=True)
    #pyplot.plot(x, vac_count, save='defect_vac_diffusion_count')
    #pyplot.plot(x, int_count, save='defect_int_diffusion_count')
    return DV, DW, [substitutions, vacancies, intersticial, x]
    #find the vacancies
#if __name__ == '__main__':
#    vac_near()
if __name__ == '__main__':
    #run_debug() #run_all()
    import MD
    from MD.dna_scripts.square import drift_remove_all
    print 'MD.L is '
    M=readxyz.ReadCord(trajectory = 'surface.xyz',frames=777)
    last = 777
    MD.L = [103.85,103.85,103.84]
    MD.last = last
    try:
        MD.V = util.pickle_load('Vs.pkl')
        MD.W = util.pickle_load('Ws.pkl')
        MD.VW = util.pickle_load('VWs.pkl')
    except:
        MD.V = M.cord_auto(['V'])
        MD.W = M.cord_auto(['W'])
        MD.VW = M.cord_auto(['V','W'])
        util.pickle_dump(MD.V,'Vs.pkl')
        util.pickle_dump(MD.W,'Ws.pkl')
        util.pickle_dump(MD.VW,'VWs.pkl')
    CV = util.pickle_load('CV.pkl')
    CW = util.pickle_load('CW.pkl')
    try:
        print 'attempting to load DW, DV'
        DV = util.pickle_load('DV_s.pkl')
        DW = util.pickle_load('DW_s.pkl')
        plot_def = util.pickle_load('plot_def_s.pkl')
    except:
        print 'finding the defects'
        DV, DW, plot_def =  find_defects_mod(CW,CV,MD.VW,MD.V,MD.W,MD.L,n_finish=MD.last-1,n_start=0,delta=1)
    plots = util.pickle_load('plot_def.pkl')
    for index,k in enumerate(x):
        Vn = []
        Wn = []
        for i in V[k]:
            if abs(i[0]) < L[0]:
                Vn.append(i)
        for i in W[k]:
            if abs(i[0]) < L[0]:
                Wn.append(i)
        particle_frame.append(len(Vn)+len(Wn))
    pyplot.plot3(plots[3], plots[0], plots[3], particle_frame, plots[3],
            plots[2], label1='substitutions',label2='number of particles',label3='intersticial',
            save='defects_per_frame',showleg=True)
