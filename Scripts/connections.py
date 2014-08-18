# -*- coding: utf-8 -*-
import sys
import os
sys.path.append('/home/cdknorow/Dropbox/Software/')
sys.path.append('/home/cdknorow/Dropbox/Software/MD')
import numpy as np
import math
import matplotlib.pyplot as plt
import MD
import MD.util as util
import MD.base.points as points
import MD.plot.pyplot as pyplot
reload(MD)
reload(pyplot)
#################################################
# What to run
#################################################
def run_single():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    index = M.get_index(['C','G'])
    total = len(index)
    fid = open('conn.dat','r')
    con = fid.readlines()
    fid.close()
    fid = open('conn_fh.dat','w')
    for i,line in enumerate(con):
        fid.write(('%i %.4f\n'%(i, float(line)/total))) 
    fid.close()
    fid = open('conn_average.data','w')
    #find the average connections every delata timesteps
    av = 0
    delta =20
    count = 0
    for i,line in enumerate(con):
        av += float(line)/total
        count +=1
        if count == delta:
            fid.write(('%i %.4f\n'%(i-5, av/delta))) 
            count = 0
            av = 0
    fid.close()
def find_arms():
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    #initial file
    dirname = os.getcwd().partition('/')[-1]
    print "Starting:",dirname.split('/')[-1]
    M=MD.ReadCord()
    index_M = M.get_index(['M'])
    index_N = M.get_index(['N'])
    index_V = M.get_index(['V','W'])
    total_N = len(index_N)
    total_M = len(index_M)
    total_V = float(len(index_V))
    print 'total number of particles'
    print total_V
    print 'number of arms'
    print total_M/total_V
    print 'number of beads per NC'
    print total_N/total_V
def find_box():
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    L_cont = M.box_volume()
    L_cont[-1] = L_cont[-2]
    fid = open('box_length.txt','w')
    fid.write('frame  L\n') 
    for i in range(0,len(L_cont)):
        fid.write(('%i   %.2f\n')%(i,L_cont[i][0]))
    fid.close()

#For single directories
if __name__ == '__main__':
    find_box()
    find_arms()
    run_single()


