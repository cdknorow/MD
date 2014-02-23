import numpy as np
import MD.util as util
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
#####################################################
# Look at defects 
#####################################################
def find_lattice(MD,VW,L,n_finish=1,n_start=0,delta=20,filters=0.2):
    # Identify the bcc lattice as a grid and find the recipricol lattice
    #Find Structure Factor
    #print VW.shape
    redo = True
    while redo == True:
        stmp,qx = sf.sfactor(VW[-10:-1],L=L,l=10)
        S,Q,primitive_vectors = sf.sort_sfpeaks(stmp,qx)
        print np.array(Q)
        Q, S, primitive_vectors = sf.sf_filter(Q, S, primitive_vectors,
                filters=filters, save = 'sffilter')
        a,b = MD.analysis.sfactor.recipricol(Q,primitive_vectors)
        #Find the best fit for the qlattice at a particular frame by looking
        #for the lowest value of the msdr
        x = np.arange(n_start,n_finish,delta)
        defects = np.zeros((len(x),15))
        lowest=10000
        for k in x:
            for i in range(10):
                crystal, c = make.make_bcc(a[0], a[1], a[2], VW[k][i], L[0]+5)
                msdr_x =  msdr.msd_r(VW[k], crystal, L)
                if lowest > msdr_x:
                    lowest = msdr_x
                    index = i
                    frame = k
        crystal, cword = make.make_bcc(a[0], a[1], a[2], VW[frame][index], L[0]+5)
        print np.array(b)
        print np.array(a)
        if raw_input('Do you want to retry?[no]/yes  ')  == 'yes':
            redo = True
        else:
            redo = False
    util.pickle_dump(cword,'cword.pkl')
    util.pickle_dump(crystal,'ccrystal.pkl')
    #filter out the comensurate crystal
    C = []
    for i in crystal:
        #make sure there aren't two that are close due to pbc
        if rcut_neighbors_point(np.array(C),i,L,rcut=6) == 0:
            C.append(i)
    #write out the lattice
    fid = open('testlattice.xyz','w')
    fid.write(('%i\n')%(len(C)))
    fid.write('Atoms\n')
    for i in range(len(C)):
        fid.write(('%c   %f   %f   %f\n')%('A',C[i][0],C[i][1], C[i][2]))
    fid.close()
    util.pickle_dump(C,'C.pkl')

    return crystal, C
#filter the 'qlattice.xyz' file so that we only have the right number of atoms
#
# returns: CV CW : The real space lattice projection of the crystal
# for V and W spots.
#
#write projection of crystal
#If this isn't working check the rcut value in find_neighbors
def filter_lattice(C, Vshape, Wshape, L):
    import random
    #### find the neighbors to a point
    def find_neighbors(M,point,L,rcut=17):
        neighbors = []
        for i in range(len(M)):
            if points.dist(M[i],point,L)[0]>1:
                if points.dist(M[i],point,L)[0]<rcut:
                    neighbors.append(M[i])
        return neighbors

    CB = []
    CA = []

    point = C[0]

    CA.append([point[0],point[1],point[2]])

    for i in range(len(C)):
        neighbors = find_neighbors(C,point,L)
        point = neighbors[random.randint(0,len(neighbors)-1)]
        for i in neighbors:
            CB.append([i[0],i[1],i[2]])
        neighbors = find_neighbors(C,point,L)
        point = neighbors[random.randint(0,len(neighbors)-1)]
        for i in neighbors:
            CA.append([i[0],i[1],i[2]])

    CV = np.zeros((len(C)/2,3))
    CW = np.zeros((len(C)/2,3))
    CV += L[0]*2
    CW += L[0]*2
    print 'print info'
    print CA[0]
    print CV.shape
    print len(C)
    print len(CA)
    print Vshape

    count = 0
    for i, point in enumerate(CA):
        if rcut_neighbors_point(CV,point,L,rcut=6) == 0:
            print count
            CV[count] = point
            count += 1
    count = 0
    for i, point in enumerate(CB):
        if rcut_neighbors_point(CW,point,L,rcut=6) == 0:
            CW[count] = point
            count += 1

    print CW.shape
    print CV.shape
    #write out hte lattice
    fid = open('projlattice.xyz','w')
    fid.write(('%i\n')%(len(C)))
    fid.write('Atoms\n')
    for i in range(len(CV)):
        if CV[i][0] != 0:
            fid.write(('%c   %f   %f   %f\n')%('A',CV[i][0],CV[i][1], CV[i][2]))
    for i in range(len(CW)):
        #there is a problem with writing one to many when putting in
        #interstitials
        if CW[i][0] != 0:
            fid.write(('%c   %f   %f   %f\n')%('B',CW[i][0],CW[i][1], CW[i][2]))
    fid.close()
    util.pickle_dump(CV,'CV.pkl')
    util.pickle_dump(CW,'CW.pkl')
    return CV, CW
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
        for i in range(CV.shape[0]):
            N =  ws_neighbors_point(V[k],CV[i],L,i,rmin=rmin,rmax=rmax)[0]
            N_V.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(W[k],CV[i],L,i,rmin=rmin,rmax=rmax)[0]
                N_W.extend(N2)
                num = -len(N2)
            DV[index][i] = num
        ###########################
        for i in range(CW.shape[0]):
            N =  ws_neighbors_point(W[k],CW[i],L,i+W.shape[1],rmin=rmin,rmax=rmax)[0]
            N_W.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(V[k],CW[i],L,i+V.shape[1],rmin=rmin,rmax=rmax)[0]
                N_V.extend(N2)
                num = -len(N2)
            DW[index][i] = num
        #find the atoms that haven't been placed on the lattice yet
        IV = points.difference(range(V.shape[1]),N_V)
        IW = points.difference(range(W.shape[1]),N_W)
        ##
        print 'atoms not added listed after first sorting'
        print IV, IW
        print DW[index]
        print DV[index]
        for i in IV:
            #find closest lattice point
            nw, dw = close_neighbors_point(CW,V[k][i],L)
            nv, dv = close_neighbors_point(CV,V[k][i],L)
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
            nw, dw = close_neighbors_point(CW,W[k][i],L)
            nv, dv = close_neighbors_point(CV,W[k][i],L)
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
        IV = points.difference(range(V.shape[1]),N_V)
        IW = points.difference(range(W.shape[1]),N_W)
        ##
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
                    print pr
                    fid.write(pr)
                if A[index][i] == -1:
                    pr = 'substitution ' + str(i+C)+ '\n'
                    #print pr
                    fid.write(pr)
                if A[index][i] == -2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
                    print pr
                    fid.write(pr)
                if A[index][i] == 2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
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
    print substitutions
    print vacancies
    print intersticial
    util.pickle_dump(DV,'DV.pkl')
    util.pickle_dump(DW,'DW.pkl')
    util.pickle_dump([substitutions, vacancies, intersticial, x],'plot_def.pkl')
    pyplot.plot3(x, substitutions, x, vacancies, x, intersticial,
            label1='substitutions',label2='vacancies',label3='intersticial',
            save='defects_per_frame',showleg=True)
    pyplot.plot(x, vac_count, save='defect_vac_diffusion_count')
    pyplot.plot(x, int_count, save='defect_int_diffusion_count')
    return DV, DW, [substitutions, vacancies, intersticials, x]
    #find the vacancies
#for normal systems
def find_defects(CV,CW,VW,V,W,L,n_finish=1,n_start=0,delta=20,filters=0.2):
    from MD.analysis.nearest_neighbor import ws_neighbors_point
    from MD.analysis.nearest_neighbor import close_neighbors_point
    #######
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
    s = max(W.shape[1],V.shape[1])
    #for vacancy
    #s = max(W.shape[1],V.shape[1])
    DW = np.zeros((len(x),s))
    DV = np.zeros((len(x),s))
    vac_count = np.zeros((len(x),1))
    int_count = np.zeros((len(x),1))
    fid = open('defects.txt','w')
    master = [[],[]]
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
        #correct spot  = 1
        #vacancies = 0
        #interstitials = 2
        #substitions are negative
        for i in range(V.shape[1]):
            N =  ws_neighbors_point(V[k],CV[i],L,i,rmin=rmin,rmax=rmax)[0]
            N_V.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(W[k],CV[i],L,i,rmin=rmin,rmax=rmax)[0]
                N_W.extend(N2)
                num = -len(N2)
            DV[index][i] = num
        ###########################
        for i in range(V.shape[1]):
            N =  ws_neighbors_point(W[k],CW[i],L,i+W.shape[1],rmin=rmin,rmax=rmax)[0]
            N_W.extend(N)
            num = len(N)
            if num == 0:
                N2 = ws_neighbors_point(V[k],CW[i],L,i+W.shape[1],rmin=rmin,rmax=rmax)[0]
                N_V.extend(N2)
                num = -len(N2)
            DW[index][i] = num
        #find the atoms that haven't been placed on the lattice yet
        IV = points.difference(range(V.shape[1]),N_V)
        IW = points.difference(range(W.shape[1]),N_W)
        print 'atoms not added list for debugging'
        print IV, IW
        for i in IV:
            #find closest lattice point
            nw, dw = close_neighbors_point(CW,V[k][i],L)
            nv, dv = close_neighbors_point(CV,V[k][i],L)
            if dw <= dv:
                #check to see if there is already an atom at that point
                print index
                print nw
                print DW.shape
                print CW.shape
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
            nw, dw = close_neighbors_point(CW,W[k][i],L)
            nv, dv = close_neighbors_point(CV,W[k][i],L)
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
        IV = points.difference(range(V.shape[1]),N_V)
        IW = points.difference(range(W.shape[1]),N_W)
        ##
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
                    print pr
                    fid.write(pr)
                if A[index][i] == -1:
                    pr = 'substitution ' + str(i+C)+ '\n'
                    print pr
                    fid.write(pr)
                if A[index][i] == -2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
                    print pr
                    fid.write(pr)
                if A[index][i] == 2:
                    pr = 'interstitial ' + str(i + C)+ '\n'
                    try:
                        def_list[1].extend(i+C)
                    except:
                        def_list[1].append(i+C)
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
    print substitutions
    print vacancies
    print intersticial
    util.pickle_dump(DW,'DW.pkl')
    util.pickle_dump(DV,'DV.pkl')
    util.pickle_dump([substitutions, vacancies, intersticial,
        x],'plot_def.pkl')
    pyplot.plot3(x, substitutions, x, vacancies, x, intersticial,
            label1='substitutions',label2='vacancies',label3='intersticial',
            save='defects_per_frame',showleg=True)
    pyplot.plot(x, vac_count, save='defect_vac_diffusion_count')
    pyplot.plot(x, int_count, save='defect_int_diffusion_count')
    return DV, DW, [substitutions, vacancies, intersticial, x]
    return DV, DW
    #find the vacancies
#write the vacancey into an xyz file
def write_vac(DV,DW,CV,CW,delta=5):
    def find_vac(A,check):
        vac = []
        for i,j in enumerate(A):
            if j == check:
                vac.append(i)
        return vac
    vacA = []
    vacB = []
    for k in range(DW.shape[0]):
        vacA.append(find_vac(DV[k],0))
        vacB.append(find_vac(DW[k],0))

    fid = open('defect_vac.xyz','w')
    for k in range(DW.shape[0]):
        for i in range(delta):
            fid.write(('%i\n')%(2*DV.shape[1]))
            fid.write('Atoms\n')
            print DV.shape
            print CV.shape
            for i in range(DV.shape[1]):
                if i in vacA[k]:
                    fid.write(('O %.2f %.2f %.2f\n')%(CV[i][0],CV[i][1],CV[i][2]))
                else:
                    fid.write('O 33 33 33\n')

            for i in range(DW.shape[1]):
                if i in vacB[k]:
                    fid.write(('O %.2f %.2f %.2f\n')%(CW[i][0],CW[i][1],CW[i][2]))
                else:
                    fid.write('O 33 33 33\n')
    fid.close()
def write_vac_new(DV,DW,CV,CW,delta=5):
    def find_vac(A,check):
        vac = []
        for i,j in enumerate(A):
            if j == check:
                vac.append(i)
        return vac
    vacA = []
    vacB = []
    for k in range(DW.shape[0]):
        vacA.append(find_vac(DV[k],0))
        vacB.append(find_vac(DW[k],0))
    fid = open('defect_vac.xyz','w')
    for k in range(DW.shape[0]):
        count = 0
        for i in range(delta):
            fid.write(('%i\n')%(1))
            fid.write('Atoms\n')
            print DV.shape
            print CV.shape
            for i in range(DV.shape[1]):
                if count == 0:
                    if i in vacA[k]:
                        fid.write(('O %.2f %.2f %.2f\n')%(CV[i][0],CV[i][1],CV[i][2]))
                        count+=1
            for i in range(DW.shape[1]):
                if count == 0:
                    if i in vacB[k]:
                        fid.write(('O %.2f %.2f %.2f\n')%(CW[i][0],CW[i][1],CW[i][2]))
                        count += 1
    fid.close()
#write the vacancey into an xyz file
def write_sub(DV,DW,CV,CW,delta=1):
    def find_sub(A,check):
        sub = []
        for i,j in enumerate(A):
            if j == check:
                sub.append(i)
        return sub
    subA = []
    subB = []
    for k in range(DW.shape[0]):
        print DV[k]
        print DW[k]
        subA.append(find_sub(DV[k],1))
        subB.append(find_sub(DW[k],1))

    fid = open('defect_subs.xyz','w')
    for k in range(DW.shape[0]):
        for i in range(delta):
            fid.write(('%i\n')%(2*DV.shape[1]))
            fid.write('Atoms\n')
            for i in range(DV.shape[1]):
                if i in subA[k]:
                    fid.write(('G %.2f %.2f %.2f\n')%(CV[i][0],CV[i][1],CV[i][2]))
                else:
                    fid.write('G 100 100 100\n')

            for i in range(DW.shape[1]):
                if i in subB[k]:
                    fid.write(('R %.2f %.2f %.2f\n')%(CW[i][0],CW[i][1],CW[i][2]))
                else:
                    fid.write('R 100 100 100\n')
    fid.close()
def write_int(DV,DW,CV,CW,delta=1):
    def find_sub(A,check):
        sub = []
        for i,j in enumerate(A):
            if abs(j) == check:
                sub.append(i)
        return sub
    subA = []
    subB = []
    for k in range(DW.shape[0]):
        print DV[k]
        print DW[k]
        subA.append(find_sub(DV[k],2))
        subB.append(find_sub(DW[k],2))

    fid = open('defect_int.xyz','w')
    for k in range(DW.shape[0]):
        for i in range(delta):
            fid.write(('%i\n')%(2*DV.shape[1]))
            fid.write('Atoms\n')
            for i in range(DV.shape[1]):
                if i in subA[k]:
                    fid.write(('I %.2f %.2f %.2f\n')%(CV[i][0],CV[i][1],CV[i][2]))
                else:
                    fid.write('I 33 33 33\n')

            for i in range(DW.shape[1]):
                if i in subB[k]:
                    fid.write(('I %.2f %.2f %.2f\n')%(CW[i][0],CW[i][1],CW[i][2]))
                else:
                    fid.write('I 33 33 33\n')
    fid.close()
#find the msd of the vacancy
def vac_msd(DV,DW,CV,CW,L):
    vacancies = np.zeros((DV.shape[0],1,3))
    for k in range(len(DW)):
        print k
        #find the points that are nearest neighbor that are different
        index = np.where(DW[k] == 0)[0]
        if len(index) > 0:
            if len(index) > 1:
                index = index[1]
            vacancies[k][0]  = CW[index]
        else:
            index = np.where(DV[k] == 0)[0]
            if len(index) > 1:
                index = index[1]
            vacancies[k][0] = CV[index]
    from MD.analysis.msd import msd_no_drift as msd
    #Find the msd of the system
    x, msd = msd(vacancies, L)
    msd = msd**0.5
    util.pickle_dump([x,msd],'msd_vac.pkl')
    pyplot.plot(x,msd,xlabel='time',ylabel='msd',save='MSD_vacancey')
def int_msd(DV,DW,V,W,CV,CW,L):
    #we create an active_int which keeps track of the index at the current frame
    #we creat an all_int which keeps track of the intersticals
    #whenever the number of intersticials decerases we find the one that
    #vanished and add its trajectory to all_int
    #whenever teh number of intersticials increases we add a new trajectory 
    #to the active ints
    interstitials = np.zeros((DV.shape[0],1,3))
    last_index = 0
    active_int = []
    all_int = []
    ##########################
    # Functions
    ########################
    def find_nearest(point,CW,CV,L):
        current = 100
        for i in CW:
            d = points.dist(point,i,L)[0]
            if d < current:
                current = d
                new_point = i
        for i in CV:
            d = points.dist(point,i,L)[0]
            if d < current:
                current = d
                new_point = i
        return new_point
    #place the intersticial into the active_int array
    def place_int(new_point, active_int,list_index,L):
        current = 100
        for i in range(len(active_int)):
            if i not in list_index:
                d = points.dist(new_point,active_int[i][-1],L)[0]
                if d < current:
                    current = d
                    index = i
        try:
            index
            active_int[index].append(new_point)
            list_index.append(index)
            return index
        except:
            pass
    #find the distance between the new points and the active_int array 
    def find_int_dist(new_point, active_int,L):
        current = 100
        for i in range(len(active_int)):
            d = points.dist(new_point,active_int[i][-1],L)[0]
            if d < current:
                current = d
        return current
    def find_point(DV,CV,CW,i):
        if i < DV.shape[1]:
            new_point = CV[i]
        if i >= DV.shape[1]:
            new_point = CW[i-DV.shape[1]]
        return new_point
    ###########################################
    # Logic: Run through once
    ###########################################
    index = np.where(abs(DV[0]) == 2)[0]
    index  = np.append(index, np.where(abs(DW[0]) == 2)[0] + DW.shape[1])
    #xtrack keeps track of the frame when an interstitial is added
    x_track = []
    for i in index:
        new_point = find_point(DV,CV,CW,i)
        active_int.append([new_point])
        x_track.append(0)

        #find the point in the activt_interstitals and append it to the
        #the correct index
    last_index = index.shape[0]
    ###########################################
    #loop through the rest
    ###########################################
    for k in range(1,len(DW)):
        print k
        #find the points that are nearest neighbor that are different
        index = np.where(abs(DV[k]) == 2)[0]
        index  = np.append(index, np.where(abs(DW[k]) == 2)[0] + DW.shape[1])
        list_index = []
        #sometimes there are more than one interstitial
        #we only want to make sure we follow the same one 
        #so we look for the closest to that last position
        #if there is the same amount of intersticials
        if index.shape[0] == last_index:
            for i in index:
                new_point = find_point(DV,CV,CW,i)
                #find the point in the activt_interstitals and append it to the
                #the correct index
                place_int(new_point, active_int, list_index, L)
        #if an interstitical has been created
        if index.shape[0] > last_index:
            new_int = []
            dist_int = []
            for i in index:
                new_point = find_point(DV,CV,CW,i)
                #find the new interstitical
                dist_int.append(find_int_dist(new_point,active_int,L))
                new_int.append(new_point)
            #find the index of the max point
            max_index = dist_int.index(max(dist_int))
            for i in range(len(new_int)):
                if i != max_index:
                    place_int(new_int[i], active_int,list_index, L)
            active_int.append([new_int[max_index]])
            x_track.append(k)
        #if an intersticial has been destroyed
        if index.shape[0] < last_index:
            ann_int = []
            dist_int = []
            placed = []
            for i in index:
                new_point = find_point(DV,CV,CW,i)
                #find the annihlated interstitical
                dist_int.append(find_int_dist(new_point,active_int,L))
                ann_int.append(new_point)
            #find the index of the max point
            max_index = dist_int.index(max(dist_int))
            for i in range(len(ann_int)):
                placed.append(place_int(ann_int[i],active_int,list_index,L))
            #find the point to be removed
            dif = points.difference(range(len(active_int)),placed)
            #add the point to the list of all intersticials
            for i in dif:
                all_int.append(active_int[i])
            #remove the point from the active list
            count = 0
            for i in dif:
                del active_int[i-count]
                count+=1
        last_index = index.shape[0]
    for i in active_int:
        all_int.append(i)
    print 'len int'
    print len(all_int)
    final_int = []
    final_x = []
    #lets filter out the short lived interstitials
    for i in range(len(all_int)):
        if len(all_int[i]) > 1:
            final_int.append(all_int[i])
            final_x.append(x_track[i])
    print len(final_int)
    print 'print x'
    print len(final_x)
    from MD.analysis.msd import msd_no_drift as msd
    #Find the msd of each interstitial
    msd_all = []
    x_all = []
    for i in final_int:
        A = np.zeros((len(i),1,3))
        for j in range(len(i)):
            A[j][0][0] = i[j][0]
            A[j][0][1] = i[j][1]
            A[j][0][2] = i[j][2]
        x, msds = msd(A, L)
        msd_all.append(msds)
        x_all.append(x)
    #now we need to relate the msd to the point in time the interstital was
    #active
    msd_total = np.zeros(DV.shape[0])
    for i,k in enumerate(final_x):
        print k
        for j in range(len(msd_all[i])):
            if k+j < msd_total.shape[0]:
                msd_total[k+j] += msd_all[i][j]
    #for i in range(1,msd_total.shape[0]):
    #    msd_total[i] += msd_total[i-1]
    msd_total = msd_total**0.5
    x = range(msd_total.shape[0])
    util.pickle_dump([x,msd_total],'msd_int.pkl')
    pyplot.plot(x,msd_total,xlabel='time',ylabel='msd',save='MSD_int')
##################################333###
# identify the type atoms that make up interstitial
#########################################
# V, W are the trajecoties
# DV, DW are the defect types on lattice
# CV, CW are the lattice point locations
def int_type(V, W, DV, DW, CV, CW, L):
    #find the two nearest points in the system
    #and record their type
    def find_nearest(point, A, B, L):
        smallA = 100
        smallB = 100
        typeA = 'G'
        typeB = 'G'
        for i in range(A.shape[0]):
            d = points.dist(point,A[i],L)[0]
            if d < smallB: 
                smallA = smallB
                typeA =typeB
                smallB = d
                typeB = 'A'
            else:
                if smallA < d:
                    smallA = d
                    typeA = 'A'
        for i in range(B.shape[0]):
            d = points.dist(point,B[i],L)[0]
            if d < smallB:
                smallA = smallB
                typeA =typeB
                smallB = d
                typeB = 'B'
            else:
                if smallA < d:
                    smallA = d
                    typeA = 'B'
        return [typeA,typeB]

    totaldiff = []
    totalsame =[]
    for k in range(DV.shape[0]):
        print k
        types = []
        for i in range(DV.shape[1]):
            #if it is an interstitial, look for the 2 closests points
            if abs(DV[k][i]) == 2:
                types.append(find_nearest(CV[i],V[k],W[k],L))
        for i in range(DW.shape[1]):
            if abs(DW[k][i]) == 2:
                types.append(find_nearest(CW[i],V[k],W[k],L))
        different = 0
        same = 0
        for i in types:
            if i[0] == i[1]:
                different += 1
            else:
                same += 1
        totaldiff.append(different)
        totalsame.append(same)
    print totaldiff
    print totalsame

    x = range(V.shape[0])
    pyplot.plot2(x,totaldiff,x,totalsame,label1='different',label2='same',showleg=True,xlabel='Time',ylabel='interstistials',save='int_type')
#find what sort of defect the vacaney moves too
#find what sort of GNP moves to where the vacancy was
#count up these numbers. 
def vac_move_to(DV,DW):
    vacancies = np.zeros((DV.shape[0],1))
    #find the vacancies in the simulations
    for k in range(len(DW)):
        #find the points that are nearest neighbor that are different
        index = np.where(DW[k] == 0)[0]
        if len(index) > 0:
            if len(index) > 1:
                index = index[1]
            vacancies[k]  = index
        else:
            index = np.where(DV[k] == 0)[0]
            if len(index) > 1:
                index = index[1]
            vacancies[k] = index + DW.shape[1]
    #now that we have the indexes lets go back through and pick out where they
    #came from
    def check(vac, DW, DV):
        if vac > DW.shape[0]:
            return DV[vac-DW.shape[0]]
        else:
            return DW[vac]
    last = []
    new = []
    for k in range(1,vacancies.shape[0]):
        if vacancies[k] == vacancies[k-1]:
            last.append(0)
            new.append(0)
        else:
            last.append(check(vacancies[k][0],DW[k-1],DV[k-1]))
            new.append(check(vacancies[k-1][0],DW[k],DV[k]))
    count = 0
    count_def = 0
    print 'GNP that moved was'
    for i in last:
        if i == 1.0:
            count += 1
        if i == -1.0:
            count_def += 1
    out = open('count_def.txt','w')
    print 'not a defect'
    print count
    print 'a substitution'
    print count_def
    out.write('GNP that moved was\n')
    out.write(('not a defect %i\n')%(count))
    out.write(('a substituion %i\n')%(count_def))
    count = 0
    count_def = 0
    for i in new:
        if i == 1.0:
            count += 1
        if i == -1.0:
            count_def += 1
    print 'filled vacancy with'
    print 'not a defect'
    print count
    print 'a substitution'
    print count_def
    out.write('Filled vacancy with\n')
    out.write(('not a defect %i\n')%(count))
    out.write(('a substituion %i\n')%(count_def))
    count = 0
    count_def = 0
    print 'last vacancy filled with'
    for i in range(len(new)):
        if new[i] == 1.0:
            if new[i-1] == 0:
                count += 1
        if new[i] == -1.0:
            if new[i-1] == 0:
                count_def += 1
    print 'not a defect'
    print count
    print 'a substitution'
    print count_def
    out.write('last vacancy filled with\n')
    out.write(('not a defect %i\n')%(count))
    out.write(('a substituion %i\n')%(count_def))
    count = 0
    count_def = 0
    print 'First GNP that moved was'
    for i in range(len(last)):
        if last[i] == 1.0:
            if last[i-1] == 0:
                count += 1
        if last[i] == -1.0:
            if last[i-1] == 0:
                count_def += 1
    print 'not defect'
    print count
    print 'a substitution'
    print count_def
    out.write('First GNP moved was\n')
    out.write(('not a defect %i\n')%(count))
    out.write(('a substituion %i\n')%(count_def))

    x = range(len(new))
    #Find the msd of the system
    pyplot.plot(x,new,xlabel='Time',ylabel='vac',save='vacancey_to')
    pyplot.plot(x,last,xlabel='Time',ylabel='vac',save='vacancey_from')
def vac_surround(DV,DW,CV,CW,L):
    vacancies = np.zeros((DV.shape[0],1))
    #find the vacancies in the simulations
    for k in range(len(DW)):
        #find the points that are nearest neighbor that are different
        index = np.where(DW[k] == 0)[0]
        if len(index) > 0:
            if len(index) > 1:
                print index
                index = index[0]
            vacancies[k]  = index
        else:
            index = np.where(DV[k] == 0)[0]
            if len(index) > 1:
                index = index[0]
            vacancies[k] = index + DW.shape[1]
    #now that we have the indexes lets go back through and pick out where they
    #came from
    def check(vac, DW, DV):
        if vac > DW.shape[0]:
            return DV[vac-DW.shape[0]]
        else:
            return DW[vac]
    # find lattice point
    def find_lattice(vac, CW, CV):
        vac = int(vac[0])
        if vac > CW.shape[0]:
            return CV[vac-CW.shape[0]]
        else:
            return CW[vac]
    #find the index of the neighbors of the vacancey
    def find_neighbors(index, CW, CV, L, cut = 17):
        point = find_lattice(index,CW,CV)
        V_neighbors = []
        W_neighbors = []
        for i,j in enumerate(CV):
            if points.dist(j,point,L)[0] < cut:
                if points.dist(j,point,L)[0] > 1:
                    V_neighbors.append(i)
        for i,j in enumerate(CW):
            if points.dist(j,point,L)[0] < cut:
                if points.dist(j,point,L)[0] > 1:
                    W_neighbors.append(i)
        return V_neighbors, W_neighbors

    #Find out if the neighbors are substitutions or correctly placed
    def find_subs(V, W, DW, DV):
        subs = []
        for i in V:
            subs.append(DV[i])
        for i in W:
            subs.append(DW[i])
        return subs

    #find the number of substitutions at each step surrounding the vacancy
    subs = []
    for k in range(len(vacancies)):
        V_neigh, W_neigh = find_neighbors(vacancies[k], CW, CV, L)
        sub  = find_subs(V_neigh, W_neigh, DW[k], DV[k])
        count = 0
        for i in sub:
            if i < 0:
                count+=1
        subs.append(count)


    #find the number of substitutions surrounding a point
    #find the time that the vacancy remains there
    time = 1
    vac = []
    for k in range(2,vacancies.shape[0]):
        if vacancies[k] == vacancies[k-1]:
            time += 1
        else:
            vac.append([subs[k-1],time,k-1])
            time = 1
    out = open('vacancy_sub_time_frame','w')
    out.write('num subs, time, frame\n')
    for i in vac:
        out.write(('%i    %i    %i\n')%(i[0],i[1],i[2]))
    out.close()
    print vac
    #lets make the data accessable. How about trying a histogram
    #first sort by numer of substittuioins
    d0 = []
    d1 = []
    d2 = []
    d3 = []
    d4 = []
    for i in vac:
        if i[1] > 1:
            if i[0] == 3:
                d3.append(i[1])
            if i[0] == 2:
                d2.append(i[1])
            if i[0] == 1:
                d1.append(i[1])
            if i[0] == 0:
                d0.append(i[1])
            if i[0] == 4:
                d4.append(i[1])
    print d0
    print d1
    print d2
    print d3
    print d4

    hist0,x0=histogram_reg(d0,20)
    hist1,x1=histogram_reg(d1,20)
    hist2,x2=histogram_reg(d2,20)
    hist3,x3=histogram_reg(d3,20)
    s0 = sum(d0)/float(len(d0))
    s1 = sum(d1)/float(len(d1))
    s2 = sum(d2)/float(len(d2))
    s3 = sum(d3)/float(len(d3))

    pyplot.plot_bar([0,1,2,3],[s0,s1,s2,s3],save='vaclife')
    plt.close()

    pyplot.plot_bar(x0,hist0,save='substitutions0')
    plt.close()
    
    util.pickle_dump([d0,d1,d2,d3],'subs_lifetime.pkl')

    pyplot.plot_bar(x0,hist0,save='substitutions0')
    plt.close()
    pyplot.plot_bar(x1,hist1,save='substitutions1')
    plt.close()
    pyplot.plot_bar(x2,hist2,save='substitutions2')
    plt.close()
    pyplot.plot_bar(x3,hist3,save='substitutions3')
    plt.close()
def vac_find():
    DV = util.pickle_load('DV.pkl')
    DW = util.pickle_load('DW.pkl')
    vac_move_to(DV,DW)
def vac_near():
    import MD
    CV = util.pickle_load('CV.pkl')
    CW = util.pickle_load('CW.pkl')
    DV = util.pickle_load('DV.pkl')
    DW = util.pickle_load('DW.pkl')
    vac_surround(DV,DW,CV,CW,MD.L)
#if __name__ == '__main__':
#    vac_near()
def run_simple():
    import MD
    from MD.dna_scripts.square import drift_remove
    M=MD.ReadCord()
    L = M.box_length
    print 'MD.L is '
    print L
    last = M.frames
    L_cont = M.box_volume()
    MD.L = [L,L,L]
    MD.last = last
    try:
        MD.V = util.pickle_load('V.pkl')
    except:
        MD.V = M.cord_auto(['V'])
        MD.util.pickle_dump(MD.V,'V.pkl')
    crystal, C =  find_lattice(MD,MD.V,L_cont[MD.last-2],n_finish=MD.last-1,n_start=MD.last-2,delta=1)
def run_defect():
    #run_debug()
    #run_all()
    import MD
    from MD.dna_scripts.square import drift_remove_all
    M=MD.ReadCord()
    L = M.box_length
    print 'MD.L is '
    print L
    last = M.frames
    MD.L = [L,L,L]
    MD.last = last
    try:
        MD.V = util.pickle_load('V.pkl')
        MD.W = util.pickle_load('W.pkl')
        MD.VW = util.pickle_load('VW.pkl')
    except:
        MD.V = MD.M.cord_auto(['V'])
        MD.W = MD.M.cord_auto(['W'])
        MD.VW = MD.M.cord_auto(['V','W'])
        MD.VW, MD.V, MD.W = drift_remove_all(MD.VW[1:],MD.V[1:],MD.W[1:],MD.L)
        util.pickle_dump(MD.V,'V.pkl')
        util.pickle_dump(MD.W,'W.pkl')
        util.pickle_dump(MD.VW,'VW.pkl')
    try:
        print 'attempting to load lattice'
        crystal = util.pickle_load('ccrystal.pkl')
        cword = util.pickle_load('cword.pkl')
        C = util.pickle_load('Ccry.pkl')
    except:
        print 'finding lattice'
        crystal, C =  find_lattice(MD.VW,MD.L,n_finish=MD.last-1,n_start=MD.last-2,delta=1)
        util.pickle_dump(C,'Ccry.pkl')
        util.pickle_dump(crystal,'ccrystal.pkl')
    try:
        print 'attempting to load CV and CW'
        CV = util.pickle_load('CV.pkl')
        CW = util.pickle_load('CW.pkl')
    except:
        print 'filtering the lattice'
        CV, CW = filter_lattice(C, MD.V.shape[1], MD.W.shape[1], MD.L)
        util.pickle_dump(CV,'CW.pkl')
        util.pickle_dump(CW,'CV.pkl')
    #try:
    #    print 'attempting to load DW, DV'
    #    DV = util.pickle_load('DV.pkl')
    #    DW = util.pickle_load('DW.pkl')
    #    plot_def = util.pickle_load('plot_def.pkl')
    #except:
    #    print 'finding the defects'
    #    try:
    #        DV, DW, plot_def = find_defects(CV,CW,MD.VW,MD.V,MD.W,MD.L,n_finish=MD.last-1,n_start=0,delta=1)
    #    except:
    #        pass
    #        print 'using mod'
    #        print 'using mod'
    #        DV, DW, plot_def = find_defects_mod(CV,CW,MD.VW,MD.V,MD.W,MD.L,n_finish=MD.last-1,n_start=0,delta=1)
    ##############################################
    ########################################
    #vac_msd(DV,DW,CV,CW,MD.L)
    #int_msd(DV,DW,MD.V,MD.W,CV,CW,MD.L)
    #plt.close()
    #vac_move_to(DV,DW)
    #plt.close()
    #x = plot_def[3]*50000
    #pyplot.plot3(x, plot_def[0], plot_def[1], plot_def[2], xlabel='Time',ylabel='defects', label1='Substitutional',label2='Vacancies',label3='Intersticials',
    #        showleg=True)
    #int_type(MD.V, MD.W, DV, DW, CV, CW, MD.L)
    #write_vac_new(DV,DW,CV,CW,delta=1)
    #write_sub(DV,DW,CV,CW,delta=1)
    #write_int(DV,DW,CV,CW,delta=1)
    vac_near()
if __name__ == '__main__':
    run_simple()
    #run_defect()

#if __name__ == '__main__':
#    import MD
#    from MD.dna_scripts.square import drift_remove_all
#    #MD.V = MD.M.cord_auto(['V'])
#    #MD.W = MD.M.cord_auto(['W'])
#    #MD.VW = MD.M.cord_auto(['V','W'])
#    #MD.VW, MD.V, MD.W = drift_remove_all(MD.VW[1:],MD.V[1:],MD.W[1:],MD.L)
#    #util.pickle_dump(MD.VW,'VW.pkl')
#    print MD.L
#    V = util.pickle_load('V.pkl')
#    W = util.pickle_load('W.pkl')
#    VW = util.pickle_load('VW.pkl')
#    #crystal = util.pickle_load('ccrystal.pkl')
#    #cword = util.pickle_load('cword.pkl')
#    CV = util.pickle_load('CV.pkl')
#    CW = util.pickle_load('CW.pkl')
#    DV = util.pickle_load('DV.pkl')
#    DW = util.pickle_load('DW.pkl')
#    #CV, CW = filter_lattice(crystal, MD.V.shape[1], MD.W.shape[1], MD.L)
#    #DV, DW = find_defects_mod(CV,CW,MD.VW,MD.W,MD.V,MD.L,n_finish=MD.last-1,n_start=0,delta=1)
#    int_type(V, W, DV, DW, CV, CW, MD.L)
#    write_vac_new(DV,DW,CV,CW,delta=1)
#    write_sub(DV,DW,CV,CW,delta=1)
#    write_int(DV,DW,CV,CW,delta=1)
#if __name__ == '__main__':
#    import MD
#    from MD.dna_scripts.square import drift_remove_all
#    #MD.V = MD.M.cord_auto(['V'])
#    #MD.W = MD.M.cord_auto(['W'])
#    #MD.VW = MD.M.cord_auto(['V','W'])
#    #MD.VW, MD.V, MD.W = drift_remove_all(MD.VW[1:],MD.V[1:],MD.W[1:],MD.L)
#    #util.pickle_dump(MD.V,'V.pkl')
#    #util.pickle_dump(MD.W,'W.pkl')
#    try:
#        crystal = util.pickle_load('ccrystal.pkl')
#        cword = util.pickle_load('cword.pkl')
#        CV = util.pickle_load('CV.pkl')
#        CW = util.pickle_load('CW.pkl')
#        d = util.pickle_load('d.pkl')
#        DV = util.pickle_load('DV.pkl')
#        DW = util.pickle_load('DW.pkl')
#        V = util.pickle_load('V.pkl')
#        W = util.pickle_load('W.pkl')
#    except:
#        print 'could not find everythong'


#    int_msd(DV,DW,MD.V,MD.W,CV,CW,MD.L)
