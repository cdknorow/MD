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
#------------------------
#sorting strings numericvally
#--------------------------
import re
def atoi(text):
    return int(text) if text.isdigit() else text

############################################
## returns the average position over a few frames
############################################
def Average(VW,L,n_start,n_frames):
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
    fid = open('star.xyz','w')
    fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
    for i in range(Average.shape[0]):
        fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
    fid.close()
    return Average

############################################
## generates the initial configuration files in .xyz format
############################################
def gen_config(delta=5):
    #setup
    dirname = os.getcwd().partition('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    last = M.frames
    V=M.cord_auto(['V','W'])
    if V.shape[1] < 1:
        V=M.cord_auto(['A'])
    x = range(0,last-delta,delta)
    for i in x:
        A = Average(V,L_cont,i,delta)
        util.print_xyz(A,L[0],atoms=V.shape[1],save='config.%i.dat'%i)
    fid = open('center.xyz','w')
    fid.write('1\n\n')
    fid.write(('V %.2f %.2f %.2f'%(L[0]/2.,L[1]/2.,L[2]/2.)))
    fid.close()

def read_q6q4(f):
    fid = open(f,'r')
    M = fid.readlines()
    V = []
    q4 = []
    q6 = []
    for line in M:
        s = line.split()
        V.append([float(s[1]),float(s[2]),float(s[3])])
        q4.append(float(s[4]))
        q6.append(float(s[5]))
    return V, q4, q6

def write_qmap_xyz(VW,Q6,Q4,frame,rcut):
    # counter bcc, hcp, fcc, mrco, liqid
    counter = [0,0,0,0,0]
    fid = open('%s.q6q4.xyz'%frame,'w')
    bcc_cut = .05
    bcc = np.array([.45,.03 ])
    hcp_cut = .035
    hcp = np.array([.475,.1 ])
    fcc_cut = .075
    fcc = np.array([.5,.17 ])
    Liq_cut = .1
    Liq = np.array([.2,.04])
    C = np.array([200,200,200])
    fid.write('%i\n\n'%(5*len(VW)))
    Qx = [[]  for n in range(5)]
    Qy = [[]  for n in range(5)]
    for j in range(len(VW)):
        w = [False for n in range(5)]
        #bcc
        if points.dist_np(bcc,np.array([Q6[j],Q4[j]]))[0] < bcc_cut:
            w[0] = True
            counter[0] = counter[0] + 1
        #hcp
        elif points.dist_np(hcp,np.array([Q6[j],Q4[j]]))[0] < hcp_cut:
            w[1] = True
            counter[1] = counter[1] + 1
        #fcc
        elif points.dist_np(fcc,np.array([Q6[j],Q4[j]]))[0] < fcc_cut:
            w[2] = True
            counter[2] = counter[2] + 1
        #liq
        elif points.dist_np(Liq,np.array([Q6[j],Q4[j]]))[0] < Liq_cut:
            w[3] = True
            counter[3] = counter[3] + 1
        #Interface
        else:
            w[4] = True
            counter[4] = counter[4] + 1
        if w[0] == True:
            fid.write(('B %.4f %.4f %.4f\n'%(VW[j][0],
                VW[j][1],VW[j][2])))
            Qx[0].append(Q4[j])
            Qy[0].append(Q6[j])
        else:
            fid.write(('B %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
        if w[1] == True:
            fid.write(('H %.4f %.4f %.4f\n'%(VW[j][0],
                VW[j][1],VW[j][2])))
            Qx[0].append(Q4[j])
            Qy[0].append(Q6[j])
        else:
            fid.write(('H %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
        if w[2] == True:
            fid.write(('F %.4f %.4f %.4f\n'%(VW[j][0],
                VW[j][1],VW[j][2])))
            Qx[0].append(Q4[j])
            Qy[0].append(Q6[j])
        else:
            fid.write(('F %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
        if w[3] == True:
            fid.write(('L %.4f %.4f %.4f\n'%(VW[j][0],
                VW[j][1],VW[j][2])))
            Qx[0].append(Q4[j])
            Qy[0].append(Q6[j])
        else:
            fid.write(('L %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
        if w[4] == True:
            fid.write(('I %.4f %.4f %.4f\n'%(VW[j][0],
                VW[j][1],VW[j][2])))
            Qx[0].append(Q4[j])
            Qy[0].append(Q6[j])
        else:
            fid.write(('I %.4f %.4f %.4f\n'%(C[0],C[1],C[2])))
    xlabel = 'q4_avg'
    ylabel = 'q6_avg'
    Label = ['bcc','hcp','fcc','liq','int']
    pyplot.plot_multi(Qx,Qy,Label,
            xlabel=xlabel,ylabel=ylabel,save=('Qmap_frame%i_rcut%.2f'%(frame,rcut)),
            limx=(0,.2),limy=(0,.6),make_marker=True)
    fid.close()
    return counter

def _run(rcut, delta=5):
    home = os.getcwd()
    binput = "/home/cdknorow/Dropbox/Software/bondorder/input"
    boutput = "/home/cdknorow/Dropbox/Software/bondorder/output"
    b = "/home/cdknorow/Dropbox/Software/bondorder"
    #first remove all the files
    def remove(destination,source):
        filelist = [ f for f in os.listdir(destination) if
                f.endswith(source) ]
        for f in filelist:
            os.remove(destination+'/'+f)
    remove(binput,".dat")
    remove(boutput,".dat")
    remove(boutput,".png")
    remove(boutput,".xyz")
    #generate the config
    gen_config(delta)
    filelist = [ f for f in os.listdir(".") if
            f.endswith(".dat") ]
    import shutil
    for f in filelist:
            shutil.move(f,binput)
    os.chdir(b)
    import auto
    auto._run(rcut)
    os.chdir('../output')
    crystals = [[],[],[],[],[]]
    dirs = os.listdir(os.getcwd())
    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    dirs.sort(key=natural_keys)
    for f in dirs:
        if f[0:3] == 'con':
            VW,Q4,Q6 = read_q6q4(f)
            counter = write_qmap_xyz(VW,Q6,Q4,float(f.split('.')[1]),rcut)
            for i in range(len(crystals)):
                crystals[i].append(counter[i])
    os.system(('cat $(find ./ -name "*.q6q4.xyz" | sort -V) > q6q4_all_rcut%.2f.xyz'%rcut))
    filelist = [ f for f in os.listdir(".") if
            f.endswith(".png") ]
    for f in filelist:
            shutil.move(f,home)
    shutil.move('q6q4_all_rcut%.2f.xyz'%rcut,home)
    os.chdir(home)
    x = [range(len(crystals[0])) for i in range(len(crystals))]
    Label = ['bcc','hcp','fcc','liquid','mrco']
    fid = open('crystal_count.dat','w')
    fid.write('rcut %.2f delta %i\n'%(rcut,delta))
    for i in range(len(crystals[0])):
        fid.write(('%i %i %i %i %i %i\n'%(i, crystals[0][i], crystals[1][i],
            crystals[2][i], crystals[3][i], crystals[4][i])))
    fid.close()
    pyplot.plot_multi(x,crystals,Label,
            xlabel='time',ylabel='crystal count',save=('crystals_rcut%.2f'%(rcut)))
#input delta sort
if __name__ == "__main__":
    rcut = float(sys.argv[1])
    try:
        delta = int(sys.argv[2])
    except:
        delta = 5
    _run(rcut,delta)
