import sys
import os
import pickle
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt





#set the fonts
matplotlib.rc('font', **{'size':'30'})
plt.close()
plt.figure(figsize=(20,12), dpi=100)
counter = 0
colors = ['r','g','b','k','y','c','m','r-','g-','b-','k-','y-','c-']
############################################
def nucleate(xmax=-1, delta=5, title='', index = 0):
    #crystal = ['bcc ','hcp ','fcc ','liq ','mirco ', 'solid']
    crystal = {'L':0,'I':1,'B':2,'H':3,'F':4}
    global counter
    global colors
    try:
        fid = open('crystal_count.dat','r')
        delta = int(fid.readline().split()[3])
        fid.close()
        if os.path.isfile('q6q4_all_rcut1.40.xyz'):
            fid = open('q6q4_all_rcut1.40.xyz','r')
        if os.path.isfile('q6q4_all_rcut11.75.xyz'):
            fid = open('q6q4_all_rcut11.75.xyz','r')
        if os.path.isfile('q6q4_all_rcut8.50.xyz'):
            fid = open('q6q4_all_rcut8.50.xyz','r')
        if os.path.isfile('q6q4_all_rcut8.50.xyz'):
            fid = open('q6q4_all_rcut8.50.xyz','r')
        cut =  fid.readline()
        fid.readline()
        VW = [[] for i in range(int(cut)/5)]
        count = 0
        for i, line in enumerate(fid.readlines()):
            if line == cut or line == '\n':
                count = 0
            else:
                if line.split()[1] == '200.0000':
                    pass
                else:
                    VW[count].append(crystal[line.split()[0]])
                    count += 1
        fid.close()
    except:
        print 'no pickle'
        return 0

    print VW[0]
    x = range(0,len(VW[0])*delta,delta)
    print len(x)
    C = [[] for i in range(20)]
    delta = 1
    for j in range(0,len(VW[0])-delta,delta):
        lm = 0
        lb = 0
        lh = 0
        lf = 0
        mb = 0
        mh = 0
        mf = 0
        ml = 0
        bf = 0
        bh = 0
        bl = 0
        bm = 0
        hf = 0
        hb = 0
        hm = 0
        hl = 0
        fm = 0
        fl = 0
        fh = 0
        fb = 0
        for i in range(len(VW)):
            #liquid to mirco
            if VW[i][j]==0 and VW[i][j+1] == 1:
                lm += 1
            #liquid to bcc
            if VW[i][j]==0 and VW[i][j+1] == 2:
                lb += 1
            #liquid to fcc
            if VW[i][j]==0 and VW[i][j+1] == 4:
                lf += 1
            #liquid to hcp
            if VW[i][j]==0 and VW[i][j+1] == 3:
                lh += 1
            #mirco to bcc 
            if VW[i][j]==1 and VW[i][j+1] == 2:
                mb += 1
            #mirco to fcc
            if VW[i][j]==1 and VW[i][j+1] == 4:
                mf += 1
            #mirco to hcp
            if VW[i][j]==1 and VW[i][j+1] == 3:
                mh += 1
            #mirco to hcp
            if VW[i][j]==1 and VW[i][j+1] == 0:
                ml += 1
            #bcc to hcp
            if VW[i][j]==2 and VW[i][j+1] == 3:
                bh += 1
            #bcc to fcc
            if VW[i][j]==2 and VW[i][j+1] == 4:
                bf += 1
            #bcc to mirco
            if VW[i][j]==2 and VW[i][j+1] == 1:
                bm += 1
            #bcc to liquid
            if VW[i][j]==2 and VW[i][j+1] == 0:
                bl += 1
            #hcp to fcc
            if VW[i][j]==3 and VW[i][j+1] == 4:
                hf += 1
            #hcp to bcc
            if VW[i][j]==3 and VW[i][j+1] == 2:
                hf += 1
            #hcp to mirco
            if VW[i][j]==3 and VW[i][j+1] == 1:
                hm += 1
            #hcp to liquid
            if VW[i][j]==3 and VW[i][j+1] == 0:
                hl += 1
            #fcc to liquid
            if VW[i][j]==4 and VW[i][j+1] == 0:
                fl += 1
            #fcc to mirco
            if VW[i][j]==4 and VW[i][j+1] == 1:
                fm += 1
            #fcc to bcc
            if VW[i][j]==4 and VW[i][j+1] == 2:
                fb += 1
            #fcc to hcp
            if VW[i][j]==4 and VW[i][j+1] == 3:
                fh += 1
        N = float(len(VW))
        C[0].append(lm/N)
        C[1].append(lb/N)
        C[2].append(lf/N)
        C[3].append(lh/N)
        C[4].append(mb/N)
        C[5].append(mh/N)
        C[6].append(mf/N)
        C[7].append(ml/N)
        C[8].append(bh/N)
        C[9].append(bf/N)
        C[10].append(bm/N)
        C[11].append(bl/N)
        C[12].append(hb/N)
        C[13].append(hm/N)
        C[14].append(hf/N)
        C[15].append(hl/N)
        C[16].append(fm/N)
        C[17].append(fb/N)
        C[18].append(fh/N)
        C[19].append(fl/N)
    l = ['liquid-mirco','liquid-bcc','liquid-fcc','liquid-hcp',
          'mirco-bcc','mirco-hcp','mirco-fcc','mirco-liquid',
          'bcc-hcp','bcc-fcc','bcc-mirco','bcc-liquid',
          'hcp-bcc','hcp-mirco','hcp-fcc','hcp-liquid',
          'fcc-mirco','fcc-bcc','fcc-hcp','fcc-liquid']
    colors = {'liquid':'b','mirco':'k','bcc':'y','hcp':'g','fcc':'r'}
    p = np.array([0,1,2,3])+index*4
    for j,i in enumerate(p):
        plt.plot(np.arange(len(C[i]))*5,C[i],'%s'%colors[l[i].split('-')[1]],lw=5,
                label =l[i].split('-')[1])
    plt.ylim((0,.2))
        #for i in range(len(C)):
    #    print i ,' ', l[i],' ',C[i]

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    config = {'-f':'','-x':'-1','-r':'-1','-D':'True','-l':'single','-i':3}
    options, arg =  getopt.getopt(sys.argv[1:],'f:x:r:D:l:i:')
    config.update( dict(options))
    if config['-f'] != '':
        os.chdir(config['-f'])
        title = config['-f'].split('_')[10]
        nucleate(title = title, xmax = int(config['-r']),index = int(config['-i']))
        os.chdir('../')
    elif len(arg) > 1:
        if config['-l']=='single':
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    dirs.append(i)
            for i in arg:
                f = dirs[int(i)]
                print f
                os.chdir(f)
                if f[0] =='F':
                    title = 'rho = %s'%f.split('_')[9]
                if f[0] =='s':
                    title = 'rho = %s'%f.split('_')[10]
                #if f[0] =='F':
                #    title = 'N = %s'%f.split('_')[13]
                #if f[0] =='s':
                #    title = 'N = %s'%f.split('_')[14]
                nucleate(title = title, xmax = int(config['-r']),index = int(config['-i']))
                os.chdir('../')
        else:
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    for j in sorted(os.listdir(i)):
                        if os.path.isdir(i+'/'+j):
                            dirs.append(i+'/'+j)
            for i in arg:
                f = dirs[int(i)]
                os.chdir(f)
                r = f.split('/')[1]
                print r
                if r[0] =='F':
                    title = 'rho = %s'%r.split('_')[9]
                if r[0] =='s':
                    title = 'rho = %s'%r.split('_')[10]
                if r[0] =='F':
                    title = 'N = %s'%r.split('_')[13]
                if r[0] =='s':
                    title = 'N = %s'%r.split('_')[14]
                nucleate(title = title, xmax = int(config['-r']),index = int(config['-i']))
                os.chdir('../../')
    else:
        if config['-l']=='single':
            dirs = []
            for i in sorted(os.listdir('.')):
                if os.path.isdir(i):
                    dirs.append(i)
            for i in arg:
                f = dirs[int(i)]
                print f
                os.chdir(f)
                if f[0] =='F':
                    title = 'rho = %s'%f.split('_')[9]
                if f[0] =='s':
                    title = 'rho = %s'%f.split('_')[10]
                #if f[0] =='F':
                #    title = 'N = %s'%f.split('_')[13]
                #if f[0] =='s':
                #    title = 'N = %s'%f.split('_')[14]
                nucleate(xmax = int(config['-r']),index = int(config['-i']))
                os.chdir('../')
    plt.xlabel = 'Time'
    plt.ylabel = ''
    plt.legend(loc=2)
    plt.show()


