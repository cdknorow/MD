# -*- coding: utf-8 -*-
import sys
import os

#setup
#take one frame
#number of particles per sna
NP = 561
#critical nucleus for frame 10
list_keep = [80, 240,261, 67, 122, 299, 384, 58, 111, 130, 368, 365, 266, 186,
             128,  310, 254, 345, 356, 301, 89, 83, 303, 278, 282,
             48, 263, 158,141, 373, 47]
#289

fid = open('xyzrewrite_short.xyz','r')
out = open('nucleus.xyz','w')

cords = fid.readlines()

#loop through all of the relevant NP
out.write(('%i\n')%(NP*len(list_keep)))
out.write('\n')
for particle in list_keep:
    for j in range(NP):
        s = cords[2 + particle*NP + j].split()
        line = s[0] + ' ' + s[1] + ' ' + s[2] + ' ' + s[3] + '\n' 
        out.write(line)
    #############
    # not implemented#
    #check to see if the beads is wrapped around the box
    #use the placement of VW bead to determine where it should be placed
    ################

