## \package MD
# \brief  analysis of molecular simulations in python
#
# MD is a python toolkit to analyze molecular dynamics
# trajectories in xyz format.
#
# It allows one to read molecular dynamics trajectories and access the
# atomic coordinates through numpy arrays. This provides a flexible and
# relatively fast framework for complex analysis tasks. 
#
# Citation
# --------
# When using MD in published work, please cite
#  	Knorowski, C.
# Thanks!
#
# Dependecies:
#---------------
# python >  2.4 but not 3.0+
# numpy - high priortiy
# catdcd - converts .dcd and .mol2 files into xyz format
# matplotlib - if you want to do any ploting
# scipy - for some functions
# python-graph - http://code.google.com/p/python-graph/
#
# Getting started
# ---------------
#
# Import the package::
# 
#  >>> import MD
#
#
# Examples
# --------
#
# see dna_scripts


from util.readcord import ReadCord  
from util import util
from util import readdcd
import os

if __name__ == '__main_':
    #dirname = os.getcwd().partition('/')[-1]
    #print "Starting:",dirname.split('/')[-1]
    #var = util.get_path_variables()
    #D = readdcd.DCD('dna.dcd')
    #last = D.numframes
    #L = D.box_length
    #M = ReadCord()
    #V = M.cord_auto(['V'])
    #W = M.cord_auto(['W'])
    #VW = M.cord_auto(['V','W'])
    print 'variables  can be loaded by'
    print 'V = M.cord_auto(["V"])'
    print 'VW = M.cord_auto(["V","W"])'
    print 'box size is L'
    print 'number of frames is last'
    print 'path variables is var'



