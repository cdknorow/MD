## \package MD.analysis.cubatic_order
# \brief This scirpt implements an alogrith to find uniaxial and biaxial ordering
#of cubes as defined in ref: http://dx.doi.org/10.1063/1.1711594
import sys
import numpy as np
import MD.base.points as points
# M [cube][x,y,z]
# Z [edge][x,y,z], edge 1-6 represent first cube
def cube_vectors(V,E,L):
    vec = []
    for i in E:
        vec.append(points.unit(points.vector1d(V,i,L)))
    return vec
#kronecker delta function
def delta(i,j):
    if i == j:
        return 1
    else:
        return 0
#kronecker delta function
def Mab(Q):
    c = max(Q[1][1],Q[0][0],Q[2][2])
    if Q[0][0]==c:
        Q[0][0] -= 0.70
        Q[1][1] -= 0.15
        Q[2][2] -= 0.15
        return 1
    if Q[1][1]==c:
        Q[1][1] -= 0.70
        Q[0][0] -= 0.15
        Q[2][2] -= 0.15
        return 2
    if Q[2][2]==c:
        Q[2][2] -= 0.70
        Q[0][0] -= 0.15
        Q[1][1] -= 0.15
        return 3
# \brief find uniaxial and biaxial ordering of cubes within a crystal
#
# \return Q00 and Q22 bond ordering parameter as defined in ref 
# http://dx.doi.org/10.1063/1.1711594
#
# \param M[frame][atoms][x,y,z]: center atom of cube
# \param Z[frame][side][x,y,z]: center points for sides of cube 
# \param L length of box
# \param select reference cube for definiing intial x,y and z axis
def old_cubatic_order(M,Z,L,select=0):
    #Choose references cube 
    ref = M[0][select]
    us  = points.unit(points.vector1d(ref,Z[0][select*6],L))
    ws  = points.unit(points.vector1d(ref,Z[0][select*6+1],L))
    vs  = points.unit(points.vector1d(ref,Z[0][select*6+4],L))
    ####################################################
    #Create u,w,v based on best alignment with reference cube
    u = []
    v = []
    w = []
    for k in range(M.shape[0]):
        for i in range(M.shape[1]):
            if i != select:
                vc = cube_vectors(M[k][i],Z[k][i*6:i*6+6],L)
                usmall = 0
                vsmall = 0
                wsmall = 0
                for j in range(1,len(vc)):
                    if points.angle_between(us,vc[j])<=points.angle_between(us,vc[usmall]):
                        usmall = j
                    if points.angle_between(ws,vc[j])<=points.angle_between(ws,vc[wsmall]):
                        wsmall = j
                    if points.angle_between(vs,vc[j])<=points.angle_between(vs,vc[vsmall]):
                        vsmall = j
                u.append(vc[usmall])
                v.append(vc[vsmall])
                w.append(vc[wsmall])
            else:
                u.append(ws)
                v.append(us)
                w.append(vs)
    #Find Ordering Tesnors
    Quu = np.zeros((3,3))
    Qww = np.zeros((3,3))
    Qvv = np.zeros((3,3))
    for N in range(len(u)):
        for i in range(3):
            for j in range(3):
                Quu[i,j] += 3. * u[N][i]*u[N][j] - delta(i,j)
                Qvv[i,j] += 3. * v[N][i]*v[N][j] - delta(i,j)
                Qww[i,j] += 3. * w[N][i]*w[N][j] - delta(i,j)
    Quu = Quu / (2.*len(u))
    Qvv = Qvv / (2.*len(v))
    Qww = Qww / (2.*len(w))
    #Find eigenvalues
    e_w = np.linalg.eig(Qww)
    e_u = np.linalg.eig(Quu)
    e_v = np.linalg.eig(Qvv)
    #Identify Laboratory Z,Y and X axis
    #print e_w,'\n',e_u,'\n',e_v
    def find_max(e):
        evalue = e[0][0]
        index = 0
        if e[0][1] > evalue:
            evalue = e[0][1]
            index = 1
        if e[0][2] > evalue:
            evalue = e[0][2]
            index = 2
        return evalue, index
    #get index of eigenvector and evalue
    e_plus_v = find_max(e_v)
    e_plus_u = find_max(e_u)
    e_plus_w = find_max(e_w)
    #find Z
    s = []
    if e_plus_v[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_v[1][e_plus_v[1]]
        Qz = Qvv
        s.append(0)
    if e_plus_v[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_v[1][e_plus_v[1]]
        Qx = Qvv
        s.append(0)
    if e_plus_u[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_u[1][e_plus_u[1]]
        Qz = Quu
        s.append(1)
    if e_plus_u[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_u[1][e_plus_u[1]]
        Qx = Quu
        s.append(1)
    if e_plus_w[0] == max(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        z = e_w[1][e_plus_w[1]]
        Qz = Qww
        s.append(2)
    if e_plus_w[0] == min(e_plus_v[0],e_plus_u[0],e_plus_w[0]):
        x = e_w[1][e_plus_w[1]]
        Qx = Qww
        s.append(2)
    if len(s) != 2:
        print 'error too Qx and Qz were not selected properly'
        asdfasdf

    if 0 not in s:
        y = e_v[1][e_plus_v[1]]
        Qy = Qvv
        s.append(0)
    if 1 not in s:
        y = e_u[1][e_plus_u[1]]
        Qy = Quu
        s.append(1)
    if 2 not in s:
        y = e_w[1][e_plus_w[1]]
        Qy = Qww
        s.append(2)

    x = np.cross(z,y)
    #print 'angles between Laboratory'
    print points.angle_between(z,y)

    #Find Q22 and Q00
    Q00 = np.dot(np.dot(z,Qz),z)

    #Q00 = (np.dot(np.dot(z,Qz),z)+np.dot(np.dot(y,Qy),y)+np.dot(np.dot(x,Qx),x))/3.
    Q22 = (np.dot(np.dot(x,Qx),x) + np.dot(np.dot(y,Qy),y) -
            np.dot(np.dot(x,Qy),x) - np.dot(np.dot(y,Qx),y))/3.

    return Q00, Q22
# \brief find uniaxial and biaxial ordering of cubes within a crystal
#
# \return Q00 and Q22 bond ordering parameter as defined in ref 
# http://dx.doi.org/10.1063/1.1711594
#
# \param M[frame][atoms][x,y,z]: center atom of cube
# \param Z[frame][side][x,y,z]: center points for sides of cube 
# \param L length of box
# \param select reference cube for definiing intial x,y and z axis
def cubatic_order(M,Z,L,frame=1):
    #us = np.array(([1/3.**0.5,1/3.**0.5,1/3.**0.5]))
    #ws = np.array(([1/2.**0.5,0,-1/2.**0.5]))
    #vs = np.array(([1/6.**0.5,-2/6.**0.5,1/6.**0.5]))
    #us = np.array(([1,0,0]))
    #ws = np.array(([0,1,0]))
    #vs = np.array(([0,0,1]))
    u_vectors = []
    u_vectors.append(np.array(([1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
    u_vectors.append(np.array(([-1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
    u_vectors.append(np.array(([0, 1/2.**0.5, 1/2.**0.5])))
    u_vectors.append(np.array(([1/2.**0.5, 0, 1/2.**0.5])))
    u_vectors.append(np.array(([1/2.**0.5, 1/2.**0.5, 0])))
    u_vectors.append(np.array(([-1/2.**0.5, 1/2.**0.5, 0])))
    u_vectors.append(np.array(([-1/2.**0.5, 0, 1/2.**0.5])))
    u_vectors.append(np.array(([1, 0, 0])))
    u_vectors.append(np.array(([0, 1, 0])))
    u_vectors.append(np.array(([0, 0, 1])))


    for again in range(2):
        if again == 0:
            fid = open('gaussmap.xyz','w')
            fid.write('%i\n\n'%(M.shape[0]*6))
        ####################################################
        #Create Master Vector
        U = [[] for i in range(len(u_vectors))]
        count = [[0,i] for i in range(len(u_vectors))]
        for i in range(M.shape[0]):
            vc = cube_vectors(M[i],Z[i*6:i*6+6],L)
            if again == 0:
                for i in range(len(vc)):
                    fid.write(('N %.2f %.2f %.2f\n')%(vc[i][0],vc[i][1],vc[i][2]))
            u_small = [0 for j in range(len(u_vectors))]
            for j in range(len(u_small)):
                for k in range(len(vc)):
                    theta = points.angle_between(u_vectors[j],vc[u_small[j]])
                    theta_new = points.angle_between(u_vectors[j],vc[k])
                    if theta_new<=theta:
                        u_small[j] = k 
                        t = points.angle_between(u_vectors[j],vc[u_small[j]])
                if t < 1.5:
                    U[j].append(vc[u_small[j]])
                else:
                    count[j][0]+=1
        fid.close()


        def Qab(u):
            Quu = np.zeros((3,3))
            for N in range(len(u)):
                for i in range(3):
                    for j in range(3):
                        Quu[i,j] +=  u[N][i]*u[N][j]
            Quu = Quu / (len(u))
            e_u = np.linalg.eig(Quu)
            return Quu,e_u
        eigen_values = []
        for vv,u in enumerate(U):
            Quu, e_u = Qab(u)
            Mab(Quu)
            print '----------------------'
            print 'allign vector'
            print u_vectors[vv]
            print 'Quu_ab'
            print Quu
            print 'eigan values of <n_a,n_b>'
            print e_u[0]
            #print 'eigan values of <n_a,n_b>'
            #print e_u[1]
            #print '----------------------'
            eigen_values.append(max(e_u[0]))
        #count.sort()
        #print count
        ### TO DO
        if again <1:
            m = eigen_values.index(max(eigen_values))
            new_vector = points.line_3D_fit(np.array(U[m]))
            u_vectors.append(new_vector)
    fid = open('vectors%i.xyz'%frame,'w')
    m = eigen_values.index(max(eigen_values))
    fid.write('130\n\n')
    fid.write('Z 0 0 0\n')
    fid.write(('Z %.4f %.4f %.4f\n')%(u_vectors[m][0],u_vectors[m][1],u_vectors[m][2]))
    for v in U[m]:
        fid.write(('N %.4f %.4f %.4f\n')%(v[0],v[1],v[2]))
    fid.close()
    print max(eigen_values)
    print '--'
    if max(eigen_values) == eigen_values[-1]:
        print 'new is better!'
    return (max(eigen_values)-0.7)/.3

# \brief find uniaxial and biaxial ordering of cubes within a crystal
#
# \return Q00 and Q22 bond ordering parameter as defined in ref 
# http://dx.doi.org/10.1063/1.1711594
#
# \param M[frame][atoms][x,y,z]: center atom of cube
# \param Z[frame][side][x,y,z]: center points for sides of cube 
# \param L length of box
# \param select reference cube for definiing intial x,y and z axis
def cubatic_order_eigenvalues(M,Z,L,frame=1,delta=5):
    #us = np.array(([1/3.**0.5,1/3.**0.5,1/3.**0.5]))
    #ws = np.array(([1/2.**0.5,0,-1/2.**0.5]))
    #vs = np.array(([1/6.**0.5,-2/6.**0.5,1/6.**0.5]))
    #us = np.array(([1,0,0]))
    #ws = np.array(([0,1,0]))
    #vs = np.array(([0,0,1]))
    u_vectors = []
    import math
    #north south and east poles
    #for i in range(5):
    #    u_vectors.append(np.array(([0, math.sin(i*math.pi/(6*4)),
    #        math.cos(i*math.pi/(6.*4))])))
    #u_vectors.append(np.array(([0.72, .37, .59])))
    #v = np.array([-6.86,-6.46,7.43])
    #v = np.array([44.7-56.0,2.76,25.469])
    #u_vectors.append(v/np.dot(v,v)**0.5)
    #v = np.array(([-4.59, -0.01, 11.09]))
    #u_vectors.append(v/np.dot(v,v)**0.5)
    #v = np.array(([-2.98, -2.25, 11.41]))
    #u_vectors.append(v/np.dot(v,v)**0.5)
    v = np.array(([-9.38, 7.28, -1.78]))
    u_vectors.append(v/np.dot(v,v)**0.5)
    #u_vectors.append(np.array(([1, 0, 0])))
    #u_vectors.append(np.array(([0, 1, 0])))
    #u_vectors.append(np.array(([0, 0, 1])))
    #u_vectors.append(np.array(([1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
    ##u_vectors.append(np.array(([-1/3.**0.5, 1/3.**0.5, 1/3.**0.5])))
    ###u_vectors.append(np.array(([0, 1/2.**0.5, 1/2.**0.5])))
    #u_vectors.append(np.array(([1/2.**0.5, 0, 1/2.**0.5])))
    #u_vectors.append(np.array(([1/2.**0.5, 1/2.**0.5, 0])))
    #u_vectors.append(np.array(([-1/2.**0.5, 1/2.**0.5, 0])))
    #u_vectors.append(np.array(([-1/2.**0.5, 0, 1/2.**0.5])))


    #fid = open('gaussmap%.2f.xyz'%L[0],'w')
    fid = open('gaussmap%.2f.xyz'%frame,'w')
    fid.write('%i\n\n'%(M.shape[1]*6*M.shape[0]))
    ####################################################
    #Create Master Vector
    #U = [[] for i in range(len(u_vectors))]
    #count = [[0,i] for i in range(len(u_vectors))]
    #theta = 1.4
    #for index in range(delta):
    #    for i in range(M[index].shape[0]):
    #        vc = cube_vectors(M[index][i],Z[index][i*6:i*6+6],L)
    #        for j in range(len(vc)):
    #            fid.write(('N %.2f %.2f %.2f\n')%(vc[j][0],vc[j][1],vc[j][2]))
    #            for k in range(len(u_vectors)):
    #                if points.angle_between(u_vectors[k],vc[j]) < theta:
    #                    U[k].append(vc[j])
    U = [[] for i in range(len(u_vectors))]
    count = [[0,i] for i in range(len(u_vectors))]
    for index in range(delta):
        for i in range(M[index].shape[0]):
            vc = cube_vectors(M[index][i],Z[index][i*6:i*6+6],L)
            for r in range(len(vc)):
                fid.write(('N %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            u_small = [0 for j in range(len(u_vectors))]
            for j in range(len(u_small)):
                for k in range(len(vc)):
                    theta = points.angle_between(u_vectors[j],vc[u_small[j]])
                    theta_new = points.angle_between(u_vectors[j],vc[k])
                    if theta_new<=theta:
                        u_small[j] = k 
                        t = points.angle_between(u_vectors[j],vc[u_small[j]])
                U[j].append(vc[u_small[j]])
    fid.close()
    def order_parameter(e):
        #bakos 
        b_lamda1 = (13+3*17**0.5)/36.
        b_lamda2 = 5/18.
        b_lamda3 = (13-3*17**0.5)/36.
        #isotropic
        i_lamda1 = 0.701
        i_lamda2 = 0.15
        i_lamda3 = 0.15
        e.sort()
        e = e[::-1]
        print e
        i_order = ((abs(e[0] - i_lamda1) +
                 abs(e[1] - i_lamda2) +
                 abs(e[2] - i_lamda3)) )
                 #(abs(1-i_lamda1)+i_lamda2 + i_lamda3))
        b_order = ((abs(e[0] - b_lamda1) +
                 abs(e[1] - b_lamda2) +
                 abs(e[2] - b_lamda3)) )
                #(abs(1-b_lamda1)+b_lamda2 + b_lamda3))
        #o = (abs(e[1]-i_lamda2)/(b_lamda2-i_lamda2)  +
        #    abs(e[2]-i_lamda3)/(i_lamda3-b_lamda3)  )
        return i_order/0.6

    def Nab(u):
        Quu = np.zeros((3,3))
        for N in range(len(u)):
            for i in range(3):
                for j in range(3):
                    Quu[i,j] +=  u[N][i]*u[N][j]
        Quu = Quu / (len(u))
        e_u = np.linalg.eig(Quu)
        return Quu,e_u
    order = []
    for vv,u in enumerate(U):
        Nuu, e_u = Nab(u)
        print 'vectors'
        print u_vectors[vv]
        print 'Nab'
        print Nuu
        Mab(Nuu)
        #print '----------------------'
        #print 'allign vector'
        #print u_vectors[vv]
        #print 'Quu_ab'
        #print Nuu
        #print 'eigan values of <n_a,n_b>'
        #print e_u[0]
        #print '----------------------'
        order.append(order_parameter(e_u[0]))

    #gauss map
    fid = open('vectors%i.xyz'%frame,'w')
    m = 0
    fid.write(('%i\n\n')%(len(U[0])))
    fid.write('Z 0 0 0\n')
    fid.write(('Z %.4f %.4f %.4f\n')%(u_vectors[m][0],u_vectors[m][1],u_vectors[m][2]))
    for v in U[m]:
        fid.write(('N %.4f %.4f %.4f\n')%(v[0],v[1],v[2]))
    fid.close()
    return sum(order)/len(order)

# \brief find uniaxial and biaxial ordering of cubes within a crystal
#
# \return Q00 and Q22 bond ordering parameter as defined in ref 
# http://dx.doi.org/10.1063/1.1711594
#
# \param M[frame][atoms][x,y,z]: center atom of cube
# \param Z[frame][side][x,y,z]: center points for sides of cube 
# \param L length of box
# \param select reference cube for definiing intial x,y and z axis
def cube_allign(M,Z,L):
    u_vectors = []
    #sp6 last-10
    #u_vectors.append(np.array(([-1/3**0.5, 1/3.**0.5, 1/3.**0.5])))
    #u_vectors.append(np.array(([0, 1/2**0.5, 1/2.**0.5])))
    #u_vectors.append(np.array(([-1/2**0.5, 1/2.**0.5, 0])))
    #u_vectors.append(np.array(([-1/2**0.5, 0, 1/2.**0.5])))
    #sp12 last-10
    u_vectors.append(np.array(([1/3**0.5, 1/3.**0.5, 1/3.**0.5])))
    u_vectors.append(np.array(([0, 1/2**0.5, 1/2.**0.5])))
    u_vectors.append(np.array(([1/2**0.5, 1/2.**0.5, 0])))
    u_vectors.append(np.array(([1/2**0.5, 0, 1/2.**0.5])))

    U = []
    for i in range(M.shape[0]):
        vc = cube_vectors(M[i],Z[i*6:i*6+6],L)
        u_small = 0
        for k in range(len(vc)):
            theta = points.angle_between(u_vectors[0],vc[u_small])
            theta_new = points.angle_between(u_vectors[0],vc[k])
            if theta_new<=theta:
                u_small = k
                t = points.angle_between(u_vectors[0],vc[u_small])
        U.append(vc[u_small])
    ###############################################################
    assign_vectors = []
    for i in range(len(U)):
        u_small = 0
        theta = points.angle_between(u_vectors[u_small],U[i])
        for k in range(len(u_vectors)):
            theta_new = points.angle_between(u_vectors[k],U[i])
            if theta_new<=theta:
                u_small = k
                theta = theta_new
        assign_vectors.append(u_small)
    #write to vmd an arrow to draw
    fid2 = open('color_index.txt','w')
    fid = open('arrows.tcl','w')
    fid.write('proc vmd_draw_arrow {mol start end} { \n '+
              'set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]] \n' +
               'graphics $mol cylinder $start $middle radius 0.5\n' + 
               'graphics $mol cone $middle $end radius 0.75\n}\n\n')

    def draw_arrow(fid, VW, u_vectors, index):
            x = 5 * u_vectors[0] + VW[i][0]
            y = 5 * u_vectors[1] + VW[i][1]
            z = 5 *  u_vectors[2] + VW[i][2]
            fid.write(('draw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n'
                                   %(VW[i][0], VW[i][1], VW[i][2], x, y, z)))
    for i in range(len(assign_vectors)):
            if assign_vectors[i] == 0:
                fid.write('draw color blue\n')
                fid2.write('0\n')
                draw_arrow(fid, M,u_vectors[0],i)
            if assign_vectors[i] == 1:
                fid.write('draw color red\n')
                fid2.write('1\n')
                draw_arrow(fid, M,u_vectors[1],i)
            if assign_vectors[i] == 2:
                fid.write('draw color green\n')
                fid2.write('2\n')
                draw_arrow(fid, M,u_vectors[2],i)
            if assign_vectors[i] == 3:
                fid.write('draw color orange\n')
                fid2.write('3\n')
                draw_arrow(fid, M,u_vectors[3],i)

    fid.close()

## \brief find the least fit of cubes to bokos cubes
#
#
def least_fit_cubes(M,Z,L,frame=1,bakos_pick=True):
    bakos1 =np.array([[.66, -.66, -.33],
                     [-.66, .66, .33],
                     [.66, .33, .66],
                     [-.66, -.33, -.66],
                     [-.33, -.66, .66],
                     [.33, .66, -.66]])
    bakos2 =np.array([[ -.66, .66, -.33],
                     [.66, -.66, .33],
                     [-.66, -.33, .66],
                     [.66, .33, -.66],
                     [.33, .66, .66],
                     [-.33, -.66, -.66]])
    bakos3 =np.array([[.66, .66, -.33],
                     [-.66, -.66, .33],
                     [-.33, .66, .66],
                     [.33, -.66, -.66],
                     [.66, -.33, .66],
                     [-.66, .33, -.66]])
    bakos4 = np.array([[-.66, -.66, -.33],
                     [.66, .66, .33],
                     [.33, -.66, .66],
                     [-.33, .66, -.66],
                     [-.66, .33, .66],
                     [.66, -.33, -.66]])

    ##new cubes
    #cube1 = np.array([  [-0.439 ,0.875 ,-0.206 ],
    #                    [ 0.427,-0.857 ,0.287 ],
    #                    [0.67 ,0.12 ,-0.732 ],
    #                    [-0.647 , -0.07 , 0.759 ],
    #                    [-0.597 , -0.529 , -0.603 ],
    #                    [ 0.609,0.577 ,0.544 ]])
    #cube2 = np.array([  [ -0.8590, 0.425,-0.28 ],
    #                    [0.9 ,-0.358 ,0.213 ],
    #                    [ 0.132, 0.708, 0.69],
    #                    [ 0,-.66 ,-.756 ],
    #                    [ -0.484, -0.633,0.601 ],
    #                    [ 0.489, 0.66, -0.56]])
        #rotation of bakos cubes if necessariy
    R = points.rotation_matrix(np.array([1,0,0]),0)
    fid = open('bakos.xyz','w')
    fid.write('24\n\n')
    for i in range(bakos1.shape[0]):
        bakos1[i] = np.dot(R,bakos1[i])
        fid.write('A %.3f %.3f %.3f\n'%(bakos1[i][0],bakos1[i][1],bakos1[i][2]))
    for i in range(bakos2.shape[0]):
        bakos2[i] = np.dot(R,bakos2[i])
        fid.write('B %.3f %.3f %.3f\n'%(bakos2[i][0],bakos2[i][1],bakos2[i][2]))
    for i in range(bakos3.shape[0]):
        bakos3[i] = np.dot(R,bakos3[i])
        fid.write('C %.3f %.3f %.3f\n'%(bakos3[i][0],bakos3[i][1],bakos3[i][2]))
    for i in range(bakos4.shape[0]):
        bakos4[i] = np.dot(R,bakos4[i])
        fid.write('D %.3f %.3f %.3f\n'%(bakos4[i][0],bakos4[i][1],bakos4[i][2]))
    fid.close()
    bakos_cubes = []
    if bakos_pick==True:
        print 'bakos'
        bakos_cubes = [bakos1, bakos2, bakos3, bakos4]
    #else:
    #    print 'selecting new cubes'
    #    print bakos_pick
    #    for index in bakos_pick:
    #        bakos_cubes.append(cube_vectors(M[index],Z[index*6:index*6+6],L))
    else:
        bakos_cubes = bakos_pick
    # fit the bakos cube to the vectors of the cube
    def fit_bakos(vc):
        fit = [0,0,0,0]
        #bakos cubes
        for i,bakos in enumerate(bakos_cubes):
            #vectors of cube
            for v in vc:
                theta = -1
                #find how well the cube is fit to the chosen
                #bakos cube
                for b in bakos:
                    theta_new = np.dot(v,b)
                    if theta < theta_new:
                        theta = theta_new
                fit[i] = fit[i]+theta
        if max(fit)<5.4:
            return 4
        else:
            return fit.index(max(fit))

    fid = open('animate/bakos_index%i.txt'%frame,'w')
    #fid2 = open('gauss_map_coloured%i.xyz'%frame,'w')
    #fid2.write('%i\n\n'%(M.shape[0]*6))
    #Bead = ['A','E','C','D','N']
    Bead = ['N','N','N','N','N']
    for i in range(M.shape[0]):
        index = fit_bakos(cube_vectors(M[i],Z[i*6:i*6+6],L))
        #vc = cube_vectors(M[i],Z[i*6:i*6+6],L)
        #for j in range(6):
            #fid2.write('%c %.3f %.3f %.3f\n'%(Bead[index],
            #    vc[j][0],vc[j][1],vc[j][2]))
        fid.write('%i\n'%index)
    fid.close()
## \make a gauss map
#
#
def gaussmap(M,Z,L,x,delta,scale=10):
    #fid = open('gaussmap%.2f.xyz'%L[0],'w')
    fid = open('gaussmap_video.xyz','w')
    for k in x:
        fid.write('%i\n\n'%(M.shape[1]*6*delta+1))
        #fid.write('Z %.2f %.2f %.2f\n'%(L[k][0],0.0,0.0))
        fid.write('Z %.2f %.2f %.2f\n'%(0,0.0,0.0))
        for index in range(k,k+delta):
            for i in range(M[index].shape[0]):
                vc = cube_vectors(M[index][i],Z[index][i*6:i*6+6],L[index])
                for r in range(len(vc)):
                    vc[r] =vc[r]*scale
                    #vc[r][0]+=L[index][0]
                    fid.write(('N %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
    fid.close()
## \make a gauss map
#
#
def gaussmap_color(M,Z,L,index,A,delta=1):
    #fid = open('gaussmap%.2f.xyz'%L[0],'w')
    fid = open('animate/gaussmap_color%i.xyz'%index,'w')
    fid.write('%i\n\n'%(M.shape[1]*6*+1))
    fid.write('Z %.2f %.2f %.2f\n'%(L[index][0],0.0,0.0))
    for i in range(M[index].shape[0]):
        vc = cube_vectors(M[index][i],Z[index][i*6:i*6+6],L[index])*delta
        for r in range(len(vc)):
            #vc[r] =vc[r]*L[index][0]/2.
            #vc[r][0]+=L[index][0]
            if A[i] == 0.0:
                fid.write(('B %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            if A[i] == 1.0:
                fid.write(('R %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            if A[i] == 2.0:
                fid.write(('G %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            if A[i] == 3.0:
                fid.write(('O %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            if A[i] == 4.0:
                fid.write(('P %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
    fid.close()
# \generate gaussmap with clusters
def gaussmap_cluster(M,Z,L,x,delta,scale=10):
    import MD.base.clusters as clusters
    fid = open('gaussmap_cluster.xyz','w')
    cluster_count=6
    for k in x:
        fid.write('%i\n\n'%(M.shape[1]*6*delta+1+cluster_count))
        #fid.write('Z %.2f %.2f %.2f\n'%(L[k][0],0.0,0.0))
        fid.write('Z %.2f %.2f %.2f\n'%(0,0.0,0.0))
        for index in range(k,k+delta):
            gauss_map = []
            for i in range(M[index].shape[0]):
                vc = cube_vectors(M[index][i],Z[index][i*6:i*6+6],L[index])
                for r in range(len(vc)):
                    vc[r] =vc[r]*scale
                    gauss_map.append(vc[r])
                    fid.write(('N %.2f %.2f %.2f\n')%(vc[r][0],vc[r][1],vc[r][2]))
            CM = clusters.cluster(gauss_map,rcut=2,cluster_cut=10)
            for i in CM:
                fid.write(('R %.2f %.2f %.2f\n')%(i[0],i[1],i[2]))
    fid.close()
