# -*- coding: utf-8 -*-
import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
import numpy as np
import math
import MD
#################################################
# What to run
#################################################
def write_properties():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    index_V = M.get_index(['V','W'])
    index_A = M.get_index(['A'])
    total_A = len(index_A)
    total_V = len(index_V)
    if total_A > 0:
        print 'total number of particles'
        print total_A
        total = total_A
    if total_V > 0:
        print 'total number of particles'
        print total_V
        total = total_V
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    fid = open('box_length.txt','w')
    fid.write('frame  L Nparticles\n') 
    for i in range(0,len(L_cont)):
        fid.write(('%i   %.2f    %i\n')%(i,L_cont[i][0],total))
    fid.close()

#For single directories
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
                write_properties()
            except:
                'directory failed',directory
                pass
            print '##\nfinished with:',directory

if __name__ == '__main__':
    write_properties()
