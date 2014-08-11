import struct
import os
import numpy
import string
import numpy as np
from numpy import fromstring


class DCD:
    def __init__(self,dcdfile):
        #variables
        #self.numframes
        #self.natoms
        #
        self.dcdfile = dcdfile
        self.file = open(dcdfile,'rb')
        f = self.file
        #abbreviations
        up = struct.unpack
        cs = struct.calcsize
        #Read dcd file header
        dat = up('i4c',f.read(cs('i4c')))
        self.posNSET = f.tell()
        self.NSET = up('i',f.read(cs('i')))[0]
        self.ISTART = up('i',f.read(cs('i')))[0]
        self.NSAVC = up('i',f.read(cs('i')))[0]
        f.read(cs('4i'))
        f.read(cs('i'))
        self.NAMNF = up('i',f.read(cs('i')))[0]
        self.DELTA = up('d',f.read(cs('d')))[0]
        f.read(cs('i'))
        f.read(cs('8i'))
        dat = up('i',f.read(cs('i')))[0]
        size = up('i',f.read(cs('i')))[0]
        self.remarks = []
        self.NTITLE = up('i',f.read(cs('i')))[0]
        for i in range(0,self.NTITLE):
            dat = up('80c',f.read(cs('80c')))
            self.remarks.append(
                string.strip(string.joinfields(dat,'')))
        #errors start here
        if up('i',f.read(cs('i')))[0] != size :
            print "DCD format error 4"
        if up('i',f.read(cs('i')))[0] != 4 :
            print "DCD format error 5"
        #number of atoms
        self.N = up('i',f.read(cs('i')))[0]
        if up('i',f.read(cs('i')))[0] != 4 :
            print "DCD format error 6"
        #boundary conditions
        #the end of the header
        self.pos1 = f.tell()
        up('i',f.read(cs('i')))
        self.box_length_x = up('d',f.read(cs('d')))[0]
        up('d',f.read(cs('d')))[0]
        self.box_length_y = up('d',f.read(cs('d')))[0]
        up('d',f.read(cs('d')))
        up('d',f.read(cs('d')))
        self.box_length_z = up('d',f.read(cs('d')))[0]
        up('i',f.read(cs('i')))
        self.boxhead = f.tell()-self.pos1
        #box length
        self.fixed = self.NAMNF
        #the length of the frame plus the header
        self.pos2 = self.pos1 + cs(repr(3*self.N)+'f6i')
        #the lenght of the data in the frame
        self.rlen = cs(repr(3*(self.N-self.NAMNF))+'f6i')
        #get to the end of the file
        f.seek(0,2)
        #number of frames
        self.numframes = (f.tell()-self.pos1)/(self.rlen +self.boxhead)
        #number of atoms
        self.numatoms = self.N
        self.box_V = [[0,0,0] for i in range(self.numframes)]
        print 'number of frames',self.numframes
        print 'box length x', self.box_length_x
        print 'box length y', self.box_length_y
        print 'box length z', self.box_length_z
        print 'number of atoms',self.N

    def getframe(self,fn):
        # Abbreviations
        up = struct.unpack
        cs = struct.calcsize
        f = self.file
        # Find the right point in the file
        if fn < -1*self.numframes or fn >= self.numframes :
            print 'frame out of bound'
            raise IndexError
        elif fn < 0 :
            return self.__getitem__(self.numframes + fn)
        else :
            #jump to the frame to get data from
            f.seek(self.pos1 + self.boxhead*(1+fn)+(fn)*self.rlen)
        # Read data
        size = up('i',f.read(cs('i')))[0]
        if size != cs(repr(self.N)+'f') :
            print "DCD format error 9"
        x = fromstring(f.read(cs(repr(self.N)+'f')),dtype=np.float32)
        size = up('i',f.read(cs('i')))[0]
        up('i',f.read(cs('i')))[0]
        if size != cs(repr(self.N)+'f') :
            print "DCD format error 10"
        y = fromstring(f.read(cs(repr(self.N)+'f')),dtype=np.float32)
        size = up('i',f.read(cs('i')))[0]
        up('i',f.read(cs('i')))[0]
        if size != cs(repr(self.N)+'f') :
            print "DCD format error 11"
        z = fromstring(f.read(cs(repr(self.N)+'f')),dtype=np.float32)
        size = up('i',f.read(cs('i')))[0]
        if size != cs(repr(self.N)+'f') :
            print "DCD format error 12"
        frame = np.transpose([x,y,z])
        return frame

    ###################
    # return the box size at each frame
    #################
    def getboxsize(self):
        # Abbreviations
        up = struct.unpack
        cs = struct.calcsize
        f = self.file
        # Find the right point in the file
        for fn in range(self.numframes):
            try:
                f.seek(self.pos1 + self.boxhead*(fn)+(fn)*self.rlen)
                up('i',f.read(cs('i')))
                self.box_V[fn][0] = up('d',f.read(cs('d')))[0]
                up('d',f.read(cs('d')))[0]
                self.box_V[fn][1] = up('d',f.read(cs('d')))[0]
                up('d',f.read(cs('d')))
                up('d',f.read(cs('d')))
                self.box_V[fn][2] = up('d',f.read(cs('d')))[0]
            except:
                print 'Number of Frames is Wrong'
                print "Number of Frames = ",fn
                aasdf
        return self.box_V
    def __len__(self):
        return self.numframes
    def __del__(self):
        self.file.close()
    def __repr__(self):
        return "< DCD " + self.dcdfile + " with " + repr(self.numframes) + " frames of " + repr(self.numatoms) + ">"


