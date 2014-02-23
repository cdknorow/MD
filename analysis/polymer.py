import sys
import os
import math
import MD.base.points as points
import numpy as np
class analyze_polymer():
    def __init__(self,M,L,frames,number_spacers = 3):
        print "getting cors"
        self.P = M.cord_frames(['M'],frames)
        self.S = M.cord_frames(['S'],frames)
        self.VW = M.cord_frames(['V','W'],frames)
        self.Z = M.cord_frames(['Z'],frames)
        self.frames = frames
        self.L = L
        print self.VW.shape
        print self.P.shape
        print self.S.shape
        print self.Z.shape
        self.rotation_map = []
        self.rotation_map_poly = []
        self.ndna = self.P.shape[1]/self.VW.shape[1]
        self.ndna_poly = self.S.shape[1]/self.VW.shape[1]
        self.spacers = self.S.shape[1]/self.VW.shape[1]/(self.ndna/2)
        self.edge = points.dist(self.VW[0][0],self.Z[0][1],self.L[0])[0]+0.8
        self.number_spacers = number_spacers
    def generate_rotation_map(self):
        w = np.array([[1,0,0],[0,1,0],[0,0,1]])
        ndna = self.ndna
        for count,k in enumerate(self.frames):
            location = []
            location_poly = []
            print k
            for i in range(self.VW.shape[1]):
                #we must rotate about a specific cube reference frame
                V = self.VW[count][i]
                x_r = points.unit(points.vector1d(V,self.Z[count][i*6+1],self.L[k]))
                y_r = points.unit(points.vector1d(V,self.Z[count][i*6+2],self.L[k]))
                z_r = points.unit(points.vector1d(V,self.Z[count][i*6+5],self.L[k]))
                v = np.array([x_r,y_r,z_r])
                R = points.reference_rotation(v,w)
                for j in range(ndna):
                    d = points.dist(V,self.P[count][j+i*ndna],self.L[k])[0]
                    c_r = points.unit(points.vector1d(V,self.P[count][j+i*ndna],self.L[k]))
                    location.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
                for j in range(self.spacers-1,self.ndna_poly,self.spacers):
                    for jj in range(self.number_spacers):
                        d = points.dist(V,self.S[count][j-jj+i*self.ndna_poly],self.L[k])[0]
                        c_r = points.unit(points.vector1d(V,self.S[count][j-jj+i*self.ndna_poly],self.L[k]))
                        location_poly.append(points.unit(np.dot(R,np.transpose(c_r)))*d)
            self.rotation_map.append(location)
            self.rotation_map_poly.append(location_poly)
    def get_polymer_stats(self,count):
        label=[]
        labelazi=[]
        zed = []
        e2e = []
        frac = []
        rho = []
        radial = []
        c2e = []
        c2s = []
        curl = [ ]
        zed_total = []
        theta = []
        theta_ideal = []
        k = self.rotation_map[count]
        for p in range(0,len(k),2):
            #find the azimuth 
            s = []
            for i in k[p]:
                s.append(abs(i))
            i = s.index(max(s))
            if i == 0:
                rho.append(points.dist(s,np.array([self.edge,0,0]),self.L[self.frames[count]])[0])
            if i == 1:
                rho.append(points.dist(s,np.array([0,self.edge,0]),self.L[self.frames[count]])[0])
            if i == 2:
                rho.append(points.dist(s,np.array([0,0,self.edge]),self.L[self.frames[count]])[0])
            zed.append(abs(k[p+1][i])-self.edge)
            e2e.append(points.dist(k[p],k[p+1],self.L[self.frames[count]])[0])
            c2e.append(points.dist(np.array([0,0,0]),k[p+1],self.L[self.frames[count]])[0])
            c2s.append(points.dist(np.array([0,0,0]),k[p],self.L[self.frames[count]])[0])
            curler = []
            for j in range(self.spacers):
                curler.append(abs(k[p+1][i])-abs(self.rotation_map_poly[count][3*p/2+self.number_spacers-(j+1)][i]))
            curl.append(min(curler))
            if curl[-1] < 0:
                sp = curler.index(min(curler))
                zed_total.append(zed[-1]+abs(curl[-1]))
                new_theta=True
                if i == 0:
                    radial.append((self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][1]**2+
                                  self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][2]**2)**0.5)
                if i == 1:
                    radial.append((self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][0]**2+
                                  self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][2]**2)**0.5)
                if i == 2:
                    radial.append((self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][0]**2+
                                  self.rotation_map_poly[count][3*p/2+self.number_spacers-(sp+1)][1]**2)**0.5)
            else:
                new_theta=False
                zed_total.append(zed[-1])
                if i == 0:
                    radial.append((k[p+1][1]**2+k[p+1][2]**2)**0.5)
                if i == 1:
                    radial.append((k[p+1][0]**2+k[p+1][2]**2)**0.5)
                if i == 2:
                    radial.append((k[p+1][0]**2+k[p+1][1]**2)**0.5)
            try:
                if new_theta == True:
                    tot_e2e = points.dist(k[p],self.rotation_map_poly[count][3*p/2+self.number_spacers-(1+i)],self.L[L_frame[count]])[0]
                    theta.append(math.asin(zed_total[-1]/tot_e2e))
                else:
                    theta.append(math.asin(zed_total[-1]/e2e[-1]))
            except:
                theta.append(math.pi/2)
            theta_ideal.append(math.pi-math.acos((c2s[-1]**2+e2e[-1]**2-c2e[-1]**2)/(2*e2e[-1]*c2s[-1])))
        return zed,zed_total,e2e,c2e,c2s,curl,theta,theta_ideal,rho,radial
        #plot some distance vs its radial distance
    def avg_radial(self,radi,y,xf):
        hist_frac_s = [[0.0,0] for i in range(len(xf))]
        for i in range(len(radi)):
            r = 0
            for j in range(len(xf)):
                if radi[i] > xf[j]:
                    r = j
            hist_frac_s[r][0] += y[i]
            hist_frac_s[r][1] += 1
        hist_frac =[]
        for i in hist_frac_s:
            if i[1] != 0:
                hist_frac.append(i[0]/i[1])
            else:
                hist_frac.append(i[0])
        return hist_frac
    #output to text file
    def write_variables(self,multi_x,multi_y,label,save='save.dat'):
        fid = open(save,'w')
        fid.write('#rho')
        for i in range(len(label)):
            fid.write('  '+label[i])
        fid.write('\n')
        for k in range(len(multi_x)):
            for i in range(len(multi_x[k])):
                fid.write("  %.2f"%multi_x[k][i])
                for j in range(len(multi_y)):
                    fid.write("  %.2f"%multi_y[j][i])
                fid.write('\n')
        fid.close()
    def write_with_theory(self,multi_x,multi_y,label,save='save.dat'):
        fid = open(save,'w')
        fid_theory = open(save+'theory.dat','w')
        fid.write('#rho')
        fid_theory.write('#rho')
        for i in range(len(label)):
            if i%2:
                fid_theory.write('  '+label[i])
            else:
                fid.write('  '+label[i])
        fid.write('\n')
        fid_theory.write('\n')
        for i in range(len(multi_x[1])):
            fid_theory.write("  %.2f"%multi_x[1][i])
            for j in range(1,len(multi_y),2):
                fid_theory.write("  %.2f"%multi_y[j][i])
            fid_theory.write('\n')

        for i in range(len(multi_x[0])):
            fid.write("  %.2f"%multi_x[0][i])
            for j in range(0,len(multi_y),2):
                fid.write("  %.2f"%multi_y[j][i])
            fid.write('\n')
        fid.close()
        fid_theory.close()
