import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.analysis.particle_distance as p_dist
import MD.util as util
import pickle
############ 
# Run from inside a folder
#########
reload(MD)
########################################3
def dump_xyz(V):
    out = open('simple.xyz','w')
    print V.shape
    for i in range(V.shape[0]):
        out.write(('%i\n\n')%(V.shape[1]))
        for j in range(V.shape[1]):
            out.write(('V %.2f %.2f %.2f\n')%(V[i][j][0],V[i][j][1],V[i][j][2]))
############################################
## returns the average position over a few frames
############################################
def Average(VW,L,n_start,n_frames,write=True,save='star.xyz'):
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
    if write:
        fid = open(save,'w')
        fid.write('%i\nL=%.4f\n'%(Average.shape[0],L[n_start][0]))
        for i in range(Average.shape[0]):
            fid.write(('V %.2f %.2f %.2f\n'%(Average[i][0],Average[i][1],Average[i][2])))
        fid.close()
    return Average
def ovitojs(i):
    fid = open('ovito.js','w')
    fid.write("node = load('dat"+str(i)+ ".imd');\n")
    fid.write("wait()\n")
    fid.write("cna = new CommonNeighborAnalysisModifier({ adaptiveMode" +
            ": true })\n")
    fid.write("node.applyModifier(cna)\n")
    fid.write('wait()\n')

    fid.write('save("data'+str(i)+'.dump", LAMMPSDumpExporter,'+
            '{columnMapping: ["Structure Type"] })\n')
    fid.close()

def crystal_count(C,delta):
    fid = open('cna_crystal_count.dat','w')
    fid.write('rcut 0 delta %i\n'%delta)
    for k in range(C.shape[0]):
        bcc = 0
        other = 0
        hcp = 0
        fcc = 0
        for i in range(C.shape[1]):
            if C[k][i] == 0:
                other += 1
            if C[k][i] == 1:
                fcc += 1
            if C[k][i] == 2:
                hcp += 1
            if C[k][i] == 3:
                bcc += 1
        fid.write('%i %i %i %i %i\n'%(k,bcc,fcc,hcp,other))
#################################################################
##print cordinates in xyz format
#################################################################
def read_lammpsdump(C,k,name='data.dump'):
    fid = open(name,'r')
    for i in range(9):
        fid.readline()
    for i, line in enumerate(fid.readlines()):
        print line
        C[k][i] == int(line.split()[0])
    fid.close()
    return C
#################################################
def get_crystals():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    delta = 10
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    last = M.frames
    try:
        VW = util.pickle_load('VW.pkl')
    except:
         VW=M.cord_auto(['V','W'])
         util.pickle_dump(VW,'VW.pkl')
    #######################################
    delta = 5
    x = range(0,last-2*delta,delta)
    C = np.zeros((len(x),VW.shape[1]),dtype=int)
    os.mkdir('imd_trajectory')
    os.chdir('imd_trajectory')
    for i,k in enumerate(x):
        A = Average(VW,L_cont,k,delta)
        util.print_imd(np.array([A]),L_cont,k)
        ovitojs(k)
        os.system('ovito --script ovito.js --nogui')
        name = 'data'+'%i'%(k)+'.dump'
        fid = open(name,'r')
        for j in range(9):
            fid.readline()
        for j, line in enumerate(fid.readlines()):
            C[i][j] = int(line.split()[0])
        fid.close()
    os.chdir('../')
    crystal_count(C,delta)



#For multiple directories
if __name__ == '__main__':
    y = []
    [y.append(x[0]) for x in os.walk(os.getcwd())]
    del y[0]
    for i in y:
        print i.split('/')[-1]
    for directory in y:
        if os.path.isdir(directory):
            os.chdir(directory)
            try:
                get_crystals()
            except:
                'directory failed',directory
                pass
            print '##\nfinished with:',directory

if __name__ == '__main__':
    get_crystals()
#    #run_compress()
