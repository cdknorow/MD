import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import powerexpansioneos as peos
import scipy.interpolate as interpolate
import h5py
reload(peos)

#set the fonts
matplotlib.rc('font', **{'size':'8'})
counter = 0
colors = ['r','g','b','k','y','c','m']
colors.extend(colors)
choose = range(9)
#choose = [10]
sigma = 1.0

def Smooth(M):
    #Find the max
    start = np.argmax(M[1:])-5 
    M_new = []
    for i in range(start):
        M_new.append(M[i])
    for i in range(start,len(M)-1):
        M_new.append((M[i-1]+M[i]+M[i+1])/3.)
    M_new.append(M[-1])
    return np.array(M_new)

def soft_potential(r,sigma=2.127279,p=9.218562):
    epsilon = 1
    return epsilon * ( 1 / (np.sinh((r/sigma))**p))

def Up_potential(r,sigma=2.127279,p=9.218562):
    return epsilon * ( sigma/r)**p

def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def potential_mean_field(title, xmin=0,xmax=8,delta=1):
    global counter
    global colors
    count = 0
    sigma = 1.
    #Au
    if xmax == 1:
        rmin=12
        rmax=22
    #NP
    if xmax == 4:
        rmin=6
        rmax=14
    if xmax == 8:
        rmin=7
        rmax=15
    if xmax == 12:
        rmin=8
        rmax=17
    bin_num=75
    #sigma = 1
    #rmin=.8
    #rmax=2.2
    ax1 = get_plot()
    print os.listdir('.')
    hpy = h5py.File('pmf_dist.hdf5','r')
    for K in sorted(hpy.keys(),key=float):
        if counter in choose:
            dist = hpy[K].value
            pmf = []
            # Get length and Number of paricles
            N = hpy[K].attrs['N_atoms']
            L = hpy[K].attrs['L']
            num_frames = hpy[K].attrs["n_frames"]
            print counter, 'L =',L/sigma,'N =',N, "Dataset",K
            print 'N_frames = ', num_frames
            L = L/sigma
            bins = np.linspace(rmin,rmax,bin_num)
            hist,r = np.histogram(dist/sigma,bins=bins)
            r_new = np.zeros(r.shape[0]-1)
            for i in range(1,r.shape[0]):
                r_new[i-1] = ((r[i]+r[i-1])/2)
            r = r_new
            r_linspace = np.linspace(r[0],r[-1],300)
            dr = r[1]-r[0]
            dhist = np.zeros(hist.shape)
            for i in range(hist.shape[0]):
                dhist[i] = hist[i] / (num_frames*N*(4*math.pi*dr*r[i]**2)*N/L**3)
            fhist = interpolate.interp1d(r, Smooth(dhist), kind='cubic')
            for i in r_linspace:
                try:
                    pmf.append(-math.log(fhist(i)))
                except:
                    pmf.append(-math.log(0.001))

            #smooth the histogram
            fid = open('hist_target%i.inp'%(count),'w')
            fid.write('#rmin rmax bins L\n')
            fid.write('#%.2f %2f %i %.4f\n'%(rmin,rmax,bin_num,L))
            for i in r_linspace:
                #interpolation might make negative numbers
                if fhist(i) > 0.001:
                    fid.write('%.8f %f\n'%(i,fhist(i)))
                else:
                    fid.write('%.8f 0\n'%(i))
            fid.close()

            if delta == 0:
                pmf = np.array(pmf)
                pmf = np.array(pmf)+abs(min(pmf))
                ax1.plot(r_linspace,pmf)
                min_index = np.argmin(pmf)
                print 'rmax =', np.linspace(r[0],r[-1],300)[min_index]
                start = 40
                coeff = peos.Soft_Potential(ax1,r_linspace[start:min_index],pmf[start:min_index],show=True)
                fid = open('pmf.inp','w')
                pmf_fit = []
                for i in range(r_linspace.shape[0]):
                    if i < start:
                        pmf_fit.append(soft_potential(r_linspace[i],coeff[0],coeff[1]))
                    elif i < min_index:
                        pmf_fit.append(pmf[i])
                    else:
                        pmf_fit.append(0)
                Force = -np.gradient(pmf_fit,r_linspace[1]-r_linspace[0])
                fid.write('#r V F\n')
                for i in range(r_linspace.shape[0]):
                    fid.write('%f %f %f\n'%(r_linspace[i],pmf_fit[i],Force[i]))
                fid.close()
                ax1.set_ylim((0,30))
                ax1.set_xlim((rmin,rmax))

            if delta == 1:
                ax1.plot(r,Smooth(dhist),'x'+colors[counter],label=K,ms=7)
                ax1.plot(r,dhist,'o'+colors[counter],label='actual',ms=4)
                f = interpolate.interp1d(r,Smooth(dhist),kind='cubic')
                r = np.linspace(r[0],r[-1],300)
                ax1.plot(r,f(r),'-'+colors[counter],label='interpolate',lw=2)
                ax1.set_ylim((0,4.5))
                ax1.set_xlim((rmin,rmax))
            count+=1
        counter += 1

def full_potential_mean_field(title, xmin=0,xmax=8,delta=1):
    global counter
    global colors
    count = 0
    sigma = 1.
    bin_num=400
    ax1 = get_plot()
    hpy = h5py.File('pmf_dist.hdf5','r')
    for K in sorted(hpy.keys(),key=float):
        if counter in choose:
            dist = hpy[K].value
            pmf = []
            # Get length and Number of paricles
            N = hpy[K].attrs['N_atoms']
            L = hpy[K].attrs['L']
            num_frames = hpy[K].attrs["n_frames"]
            print counter, 'L =',L/sigma,'N =',N, "Dataset",K
            print 'N_frames = ', num_frames
            L = L/sigma
            hist,r = np.histogram(dist/sigma,bins=bin_num)
            r_new = np.zeros(r.shape[0]-1)
            for i in range(1,r.shape[0]):
                r_new[i-1] = ((r[i]+r[i-1])/2)
            r = r_new
            r_linspace = np.linspace(r[0],r[-1],300)
            dr = r[1]-r[0]
            dhist = np.zeros(hist.shape)
            for i in range(hist.shape[0]):
                dhist[i] = hist[i] / (num_frames*N*(4*math.pi*dr*r[i]**2)*N/L**3)
            fhist = interpolate.interp1d(r, Smooth(dhist), kind='cubic')
        

            #smooth the histogram
            fid = open('full_hist_target%i.inp'%(count),'w')
            fid.write('#rmin rmax bins L\n')
            fid.write('#%.2f %2f %i %.4f\n'%(r[0],r[-1],bin_num,L))
            for i in r_linspace:
                #interpolation might make negative numbers
                if fhist(i) > 0.001:
                    fid.write('%.8f %f\n'%(i,fhist(i)))
                else:
                    fid.write('%.8f 0\n'%(i))
            fid.close()

            ax1.plot(r,Smooth(dhist),'x'+colors[counter],label=K,ms=7)
            f = interpolate.interp1d(r,Smooth(dhist),kind='cubic')
            r = np.linspace(r[0],r[-1],300)
            ax1.plot(r,f(r),'-'+colors[counter],label='interpolate',lw=2)
            ax1.set_ylim((0,3.5))
            ax1.set_xlim((10,80))
            count+=1
        counter += 1



if __name__ == "__main__":
    import getopt            ## fitted function ##
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    plt.close()
    config = {'-x':'-1','-r':'-1','-i':0}
    options, arg =  getopt.getopt(sys.argv[1:],'x:r:D:i:')
    config.update( dict(options))
    if len(arg) > 1:
        ax1 =get_plot()
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0]
            full_potential_mean_field(title = title,ax1 = ax1,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    else:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[0] 
            potential_mean_field(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    plt.legend(loc=1)
    #import carnahan
    #reload(carnahan)
    #reload(peos)
    #carnahan.starling(ax1)
    #peos.modifiedstarling(ax1,c1=1,c2=2,c3=3)
    plt.show()

