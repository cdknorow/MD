import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate as integrate 
import math
import powerexpansioneos as peos
reload(peos)
#set the fonts
matplotlib.rc('font', **{'size':'20'})
counter = 0
colors = ['r','g','b','k','y','c','m']


def get_plot():
    fig1 = plt.figure(figsize=(12,8), dpi=100)
    ax1 = fig1.add_subplot(111)
    return ax1

def WHAM(title,xmin=0,xmax=1,delta=1):
    global counter
    global colors
    def log(M, index):
        analyze = []
        for i in range(len(M)):
            analyze.append(float(M[i].split()[index]))
        return analyze
    def log_distance(M, L):
        def boundary(x1,x2,L):
            if x1-x2 > L/2:
                return ((x1-x2)-L)**2
            if x1-x2 < -L/2:
                return ((x1-x2)+L)**2
            return (x1-x2)**2
        analyze = []
        for i in range(len(M)):
            line = M[i].split()
            x1 = float(line[2])
            y1 = float(line[3])
            z1 = float(line[4])
            x2 = float(line[5])
            y2 = float(line[6])
            z2 = float(line[7])
            r2 = boundary(x1,x2,L[0])+boundary(y1,y2,L[1])+boundary(z1,z2,L[2])
            analyze.append(r2**.5)
        return analyze
    fid = open('box_length.txt','r')
    M = fid.readline()
    L = [float(M.split()[0]), float(M.split()[1]), float(M.split()[2])]
    fid.close()
    ###################################################################
    fid = open('mylog.log','r')
    M = fid.readlines()
    fid.close()
    try:
        potential = log(M[1:],M[0].split().index('pair_lj_energy'))
    except:
        potential = log(M[1:],M[0].split().index('pair_table_energy'))
    fid = open('scry_potential.dat','r')
    M = fid.readlines()
    fid.close()
    r0 = log(M[1:],M[0].split().index('d'))


    r_set = [r0[0]]
    Ni = []
    count = 1
    for i in range(1,len(r0)):
        if r0[i] != r0[i-1]:
            r_set.append(r0[i])
            Ni.append(count)
            count = 0
        count+=1
    Ni.append(Ni[-1])
    r0 = r_set
    ######################################################

    distance = log_distance(M[1:],L)
    #distacnce_sim = log(M[1:],M[0].split().index('r'))
    potential_harmonic = []
    ax1 = get_plot()
    print r0
    #plt.plot(distance,potential,'x')
    count = 0
    potential_avg = 0
    sigma=1
    for i in range(len(distance)):
        if distance[i]>3.4:
            potential_avg += potential[i]
            count+=1
    distance = np.array(distance)/sigma
    potential = np.array(potential)- potential_avg/count
    plt.plot(distance,potential,'x')
    #sigma, p, l_p, k, epsilon=0
    pmf_coef = peos.Soft_Potential(ax1,distance,potential,show=True)
    x = []
    for r in np.arange(1,3,.02):
        x.append(( 1 / (np.sinh((r/pmf_coef[0]))**pmf_coef[1])))
    print min(x)




    #get sigma^2



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
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            os.chdir(f)
            title = f.split('_')[-1]
            WHAM(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
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
            WHAM(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    delta = int(config['-i']))
            os.chdir('../')
    plt.legend(loc=3)
    #import carnahan
    #reload(carnahan)
    #reload(peos)
    #carnahan.starling(ax1)
    #peos.modifiedstarling(ax1,c1=1,c2=2,c3=3)
    plt.show()

