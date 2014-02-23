import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import shutil

#find the last atoms file to be written
def start_scripts():
    import glob
    S = glob.glob('./atoms.dump*')
    M = []
    for line in S:
        M.append(int(line.split('.')[3]))
    M.sort()
    shutil.copyfile('atoms.dump.%i.xml'%M[-1],'cont/start_dna.xml')
    shutil.copyfile('job_scritp','cont/')
#write the correct values of E and F in the inputxml script and place that script in the proper directory
def change_var(Start_L,Target_L,counter,fileread,filewrite,filename):
    data = fileread.read()
    fid = open(filewrite + filename,'w')
    s = data.replace('Change_Start',Start_L)
    s = s.replace('Change_Target',Target_L)
    s = s.replace('Change_counter',counter)
    fid.write(('%s')%(s))
############################################
## msd
############################################
def msd(VW,L,step=1):
    from MD.analysis.msd import msd
    #Find the msd of the system
    x,msd=msd(VW,L_cont[-1],step=step)
    #find the slope
    count = 0
    for i in range(-55,-45,1):
        p1 += msd[i]
        count+=1
    p1 = p1/count
    count = 0
    for i in range(-11,-1,1):
        p2 += msd[i]
        count+=1
    p2 = p2/count
    return abs(p1 - p2)/45

def analysis():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    print_box_volume(L_cont,delta=delta)
    last = M.frames
    V_index = M.get_index(['V'])
    W_index = M.get_index(['W'])
    VW=M.cord_auto(['V','W'])
    util.pickle_dump(VW,'VW.pkl')
    return msd(VW,L_cont,step=1)

#For single directories
if __name__ == '__main__':
    analysis()

