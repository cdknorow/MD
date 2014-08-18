import sys
import os
import pickle
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#set the fonts
matplotlib.rc('font', **{'size':'40'})
plt.close()
plt.figure(figsize=(20,20), dpi=100)
counter = 0
colors = ['.y','.g','.r','.b','.k','.c','.m']
############################################
def crystals(xmax=-1, time_scale=1, title='',xfit =-1,D=True, index = 3):
    global counter
    global colors
    crystal = ['bcc ','hcp ','fcc ','liq ','mirco ', 'solid']
    try:
        Q =  pickle.load(open('Qtotal.pkl','r'))
    except:
        print 'no pickles'
        return 0
    count = 0
    q = Q[index]
    for i in range(len(Q[index][0])):
        plt.plot(q[0][i],q[1][i],'%s'%colors[count], ms=10, label = crystal[count])
        count +=1

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
        crystals(title = title, xmax = int(config['-r']), xfit = int(config['-x']))
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
                crystals(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'], 
                        index = int(config['-i']))
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
                crystals(title = title, xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'],
                        index = int(config['-i']))
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
                crystals(title = '', xmax = int(config['-r']),
                        xfit = int(config['-x']), D = config['-D'], 
                        index = int(config['-i']))
                os.chdir('../')
    plt.legend(loc=4)
    plt.ylim(0.0,0.6)
    plt.xlim(0.0,0.2)
    plt.show()


