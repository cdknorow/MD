#######################################################################
# the purpose of this code is to compare be a module for analyzing point data. 
# emphasis being on all particles
#######################################################################

import sys
import os
import string
import copy
import numpy as np

##################
# Given a char, return an int with the index
##################
def cypher(c):
    A = ["A","B","C","D","E","F","G","H",'I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return A.index(c)
##################
# Given an int, return a char
##################
def acypher(i):
    A = ["A","B","C","D","E","F","G","H",'I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return A[i]
#######################################
#gets box size
#######################################
def get_box_size():
    """kind of janky way to get box size

    basically chose the last atoms dump file
    then parse the lx value """
    x=os.listdir(os.getcwd())
    x.sort()
    atom_file=[]
    for i in x:
        if i.count('atoms')==1:
            atom_file.append(i)
    print atom_file
    fid=open(atom_file[-1],'r')
    fid.readline()
    fid.readline()
    fid.readline()
    L=fid.readline()
    print L
    L=float(L.split()[1].split('"')[1])
    return L
#######################################
#get the steps where the box shrings
#######################################
def get_box_size_steps(L_cont):
        x = [0]
        L_last = L_cont[0][0]
        count = 0
        x_count = []
        for i,j in enumerate(L_cont):
            if j[0] != L_last:
                L_last = j[0]
                x.append(i)
                x_count.append(count)
                count=0
            count+=1
        x.append(len(L_cont))
        x_count.append(count)
        print x_count
        print len(x)
        print len(x_count)
        for i in range(len(x_count)):
            print i, x[i], x_count[i]
            print 'xxx'
        return x, x_count
#######################################
#returns the temp from mylog.log
#######################################
def get_temp():
    fid = open('mylog.log','r')
    log = fid.readlines()
    Temp=[]
    del log[0]
    for line in log:
        Temp.append(eval(line.split('\t')[2]))
    return np.array(Temp)
#######################################
#pulls the variables from the path name
##########################################
def get_path_variables():
    st = os.getcwd()
    st = st.split('/')[-1]
    st = st.split('_')
    var={}
    for i in range(0,len(st),2):
        #sometimes there are extra things added on so we just pass if there
        #are
        try:
            var[st[i]]=float(st[i+1]) if '.' in st[i+1] else int(st[i+1])
        except:
            pass
    return var
########################################
# find the number of frames in a dcd file
#####################################
def get_total_frames(name):
    x = os.popen('catdcd -num '+name+'.dcd')
    s = x.readlines()
    x.close()
    print s
    return int(s[6].split()[1])
####################################
# convert a dcd file to xyz
# if its there
# always flag writes a new xyz file even if one alrady exists
###################################
def dcd_to_xyz(always=False,first=1):
    #create the xyz file from dcd file
    try:
        last=get_total_frames('dna')
        print last
    except:
        print 'unable to get total number of frames'
        last=1
    if always:
        os.system(('catdcd -o xyz.txt -otype xyz -stype mol2'+ 
            ' -s dna.mol2 -first %i -last %i -dcd dna.dcd')%(first,last))
    else:
        if os.path.exists('xyz.txt')==False:
            os.system(('catdcd -o xyz.txt -otype xyz -stype mol2'+
                ' -s dna.mol2 -first %i -last %i -dcd dna.dcd')%(first,last))
        if os.path.exists('xyz.txt')==False:
            print "error no XYZ.txt file created"
            print "error no XYZ.txt file created"
            print "error no XYZ.txt file created"


#################################################################
##print cordinates in xyz format
#################################################################
def print_xyz(M,L,atoms=54,save='xyz.dat'):
    print M
    fid=open(save,'w')
    fid.write('%i\n'%(atoms))
    fid.write('%.7f\n'%(L))
    fid.write('%.7f\n'%(L))
    fid.write('%.7f\n'%(L))
    print len(M)
    for i in M:
        fid.write('%.6f %.6f %.6f\n'%(i[0]+L/2,i[1]+L/2,i[2]+L/2))
    fid.close()
#################################################################
##print cordinates in xyz format
#################################################################
def print_imd(M,L,k,save='dat'):
    fid = open(save+'%i'%k+'.imd','w')
    fid.write('#F A 1 0 0 3 0 0\n') 
    fid.write('#C number x y z\n')
    fid.write('#X %.7f 0 0\n'%L[k][0])
    fid.write('#Y 0 %.7f 0\n'%L[k][1])
    fid.write('#Z 0 0 %.7f 0\n'%L[k][2])
    fid.write('#E number x y z\n')
    for N,i in enumerate(M[0]):
        fid.write('%i %.6f %.6f %.6f\n'%(N,
            i[0]+L[k][0]/2,i[1]+L[k][1]/2,i[2]+L[k][2]/2))
    fid.close()
######################
# export data for matlab
######################
def matlab_out(x,y,name='data1'):
    fid = open(name+'.dat','w')
    for i in range(len(x)):
        fid.write(('%.4f %.4f\n')%(x[i],y[i]))
    fid.close()
#load and unload pickle files
def marshal_load(f):
    print 'loading pickle',f
    import marshal as pickle
    fid = open(f,'rb')
    S=pickle.load(fid)
    fid.close()
    return S
def marshal_dump(x,f):
    import marshal as pickle
    output = open(f,'wb')
    S=pickle.dump(x,output)
    output.close()
    return S
#load and unload pickle files
def pickle_load(f):
    print 'loading pickle',f
    import cPickle as pickle
    fid = open(f,'rb')
    S=pickle.load(fid)
    fid.close()
    return S
def pickle_dump(x,f):
    import cPickle as pickle
    output = open(f,'wb')
    S=pickle.dump(x,output)
    output.close()
    return S
# Write out the volume of the box
def print_box_volume(L_cont,delta=10):
    try:
        fid = open('box_length.txt','w')
        fid.write('frame  L\n') 
        for i in range(0,len(L_cont)):
            fid.write(('%i %.2f %.2f %.2f\n')%(i,L_cont[i][0],L_cont[i][1],L_cont[i][2]))
        fid.close()
        return 1
    except:
        return 0
