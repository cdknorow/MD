import os
import sys
import math
import numpy as np

#This is a scrip which will draw a cube given 3 perpendicular vectors and a center point. 


#\creates the sides of the cube
def write_tri(fid,c,p1,p2,p3):
    def unit(p1):
	return p1/np.dot(p1,p1)
    v1 = unit(p1)
    v2 = unit(p2)
    v3 = unit(p3)
    s1 = "{%.2f %.2f %.2f}"%(c[0]+p1[0],c[1]+p1[1],c[2]+p1[2])
    s2 = "{%.2f %.2f %.2f}"%(c[0]+p2[0],c[1]+p2[1],c[2]+p2[2])
    s3 = "{%.2f %.2f %.2f}"%(c[0]+p3[0],c[1]+p3[1],c[2]+p3[2])
    s4 = "{%.2f %.2f %.2f}"%(v1[0],v1[1],v1[2])
    s5 = "{%.2f %.2f %.2f}"%(v2[0],v2[1],v2[2])
    s6 = "{%.2f %.2f %.2f}"%(v3[0],v3[1],v3[2])
    fid.write(("draw trinorm %s %s %s %s %s %s\n")%(s1,s2,s3,s1,s2,s3))

#draw edges of the cube
def write_line(fid,c,p1,p2):
    s1 = "{%.2f %.2f %.2f}"%(c[0]+p1[0],c[1]+p1[1],c[2]+p1[2])
    s2 = "{%.2f %.2f %.2f}"%(c[0]+p2[0],c[1]+p2[1],c[2]+p2[2])
    fid.write(("draw line %s %s width 3\n")%(s1,s2))

#input # center point 
# v1,v2,v3 vectors from (0,0,0) to Zi
#
def draw_cube(fid,fid2,c,v1,v2,v3,color = 'green',material='Edgy'):
    fid2.write('V %.3f %.3f %.3f\n'%(c[0],c[1],c[2]))
    z = np.array([0,0,0])
    p1 = z -v1 -v2 -v3
    p2 = z -v1 -v2 +v3
    p3 = z -v1 +v2 +v3
    p4 = z -v1 +v2 -v3

    p5 = z + v1 -v2 -v3
    p6 = z + v1 -v2 +v3
    p7 = z + v1 +v2 +v3
    p8 = z + v1 +v2 -v3

    #fid2 = open('test.xyz','w')
    #fid2.write('8\n\n')
    #fid2.write((('N %.2f %.2f %.2f\n')%(p1[0],p1[1],p1[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p2[0],p2[1],p2[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p3[0],p3[1],p3[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p4[0],p4[1],p4[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p5[0],p5[1],p5[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p6[0],p6[1],p6[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p7[0],p7[1],p7[2])))
    #fid2.write((('N %.2f %.2f %.2f\n')%(p8[0],p8[1],p8[2])))
    #fid2.close()
    #print 'written'

    #fid.write('draw delete all\n draw color blue\n')
    fid.write('draw color '+color+'\ndraw material '+material+'\n')
    #top
    write_tri(fid,c,p7,p3,p6)
    write_tri(fid,c,p3,p2,p6)
    #bottom
    write_tri(fid,c,p8,p4,p5)
    write_tri(fid,c,p4,p1,p5)
    #left
    write_tri(fid,c,p6,p1,p2)
    write_tri(fid,c,p6,p5,p1)
    #right
    write_tri(fid,c,p7,p8,p4)
    write_tri(fid,c,p7,p4,p3)
    #front
    write_tri(fid,c,p6,p5,p8)
    write_tri(fid,c,p6,p8,p7)
    #back
    write_tri(fid,c,p2,p1,p4)
    write_tri(fid,c,p2,p4,p3)

    #Draw Edges so they appear darker
    fid.write('draw color white\n')
    write_line(fid,c,p1,p2)
    write_line(fid,c,p2,p3)
    write_line(fid,c,p3,p4)
    write_line(fid,c,p4,p1)

    write_line(fid,c,p5,p6)
    write_line(fid,c,p6,p7)
    write_line(fid,c,p7,p8)
    write_line(fid,c,p8,p5)

    write_line(fid,c,p1,p5)
    write_line(fid,c,p2,p6)
    write_line(fid,c,p3,p7)
    write_line(fid,c,p4,p8)

if __name__ == "__main__":
    fid = open('square.tcl','w')
    fid2 = open('unit_square.tcl','w')
    #center point of the cube
    C =  []
    #C.append(np.array([-0,-0,0]))
    #C.append(np.array([-0,-0,0]))
    #C.append(np.array([5,5,5]))
    #C.append(np.array([-5,-5,-5]))
    #C.append(np.array([-5,5,5]))
    #C.append(np.array([5,-5,5]))
    #C.append(np.array([5,5,-5]))
    #C.append(np.array([-5,-5,5]))
    #C.append(np.array([5,-5,-5]))
    #C.append(np.array([-5,5,-5]))
    c = np.array([0,0,0])
    R = 3
    #three perpendicular vectors to draw the cube
    #v1 = np.array([2.5,0,0])
    #v2 = np.array([0,2.5,0])
    #v3 = np.array([0,0,2.5])
    p1 = np.array([-7.21,-5.13,8.8])
    v1 = p1 - np.array([-5.97,-2.53,3.54])
    v2 = p1 - np.array([-1.55,-7.09,9.17])
    v3 = p1 - np.array([-5.6,-0.1,11.67])
    v1 = v1/np.dot(v1,v1)*R
    v2 = v2/np.dot(v2,v2)*R
    v3 = v3/np.dot(v3,v3)*R
    c = p1
    v1 =np.array([-0.84,0.54,-0.12])
    v2 =np.array([0.3,0.26,-0.92])
    v3 =np.array([0.46,0.8,0.38])
    v1 = v1*R
    v2 = v2*R
    v3 = v3*R
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([0.,0.,0.])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([7.21,5.13,-8.8])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([0.33,-4.04,-14.2])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([-4.43,7.52,-9.54])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([4.75,-11.55,-3.67])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([-6.89,-9.17,-4.4])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([11.64,-2.38,0.73])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([-11.64,2.38,-0.73])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([6.89,9.17,4.4])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([-4.75,11.55,3.67])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([4.43,-7.52,9.54])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    #c = np.array([-0.33,4.04,13.2])
    #draw_cube(fid,c,v1,v2,v3,color='red')
    C = np.array([0.,0.,0.])
    #r1 = np.array([-7.21,-5.13,8.80])
    #r2 = np.array([11.64,-2.38,0.73])
    #r3 = np.array([6.89,9.17,4.40])
    r1 = np.array([-7.492,-8.919,-4.281])
    r2 = np.array([-7.849,-4.995,8.563])
    r3 = np.array([-12.666,2.319,-0.714])
    c = C
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r1
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r2
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r3
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r1 +  r2
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r1 +  r2+r3
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r1+r3
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    c = r2+r3
    draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r2
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = r1 -  r2
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = r1 +  r2-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = r1-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = r2-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1 -  r2
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1 +  r2-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r2-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1 +  r2
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1 -  r2-r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1+r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r2+r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = -r1 +  r2 +r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    #c = r1 -  r2 +r3
    #draw_cube(fid,fid2,c,v1,v2,v3,color='red')
    for i in range(1):
        pass
        #three perpendicular vectors to draw the cube
        #v1 = np.array([3,0,0])
        #v2 = np.array([0,3,0])
        #v3 = np.array([0,0,3])
        #C = np.array([0,0,0])
        #bakos cube
        #c = np.array([0,0,-35])
        #v1 = np.array([.66,-.66,-.33])
        #v2 = np.array([.66,.33,.66])
        #v3 = np.array([-.33,-.66,.66])
        #v1 = v1*R
        #v2 = v2*R
        #v3 = v3*R
        #draw the cube
        #draw_cube(fid,c,v1,v2,v3,color='red')
        #c = np.array([0,-15,-35])
        #v1 = np.array([-.66,.66,-.33])
        #v2 = np.array([.66,.33,-.66])
        #v3 = np.array([.33,.66,.66])
        #v1 = v1*R
        #v2 = v2*R
        #v3 = v3*R
        #draw_cube(fid,c,v1,v2,v3,'green')
        #v1 = np.array([.66,.66,-.33])
        #v2 = np.array([-.33,.66,.66])
        #v3 = np.array([.66,-.33,.66])
        #v1 = v1*R
        #v2 = v2*R
        #v3 = v3*R
        #c = np.array([0,15,-35])
        #draw_cube(fid,c,v1,v2,v3,'blue')
        #v1 = np.array([-.66,-.66,-.33])
        #v2 = np.array([.33,-.66,.66])
        #v3 = np.array([-.66,.33,.66])
        #v1 = v1*R
        #v2 = v2*R
        #v3 = v3*R
        #c = np.array([0,-30,-35])
        #draw_cube(fid,c,v1,v2,v3,'orange')
    fid.close()
