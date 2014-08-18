import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#set the fonts
matplotlib.rc('font', **{'size':'20'})
plt.close()
plt.figure(figsize=(20,12), dpi=100)
counter = 0
colors = ['r','g','b','k','y','c','m']
pvt_avg = [[],[]]

def pvt(title,xmin=0,xmax=1,index=0):
    global counter
    global colors
    global pvt_avg
    def log(M, index, delta=5):
        analyze = [[] for i in range(len(M[0].split()))]
        avg_count = 0
        analyze_avg = [0]
        for line in M:
            for i in range(len(line.split())):
                analyze[i].append(float(line.split()[i]))
            analyze_avg[-1] = analyze_avg[-1] + analyze[index][-1]/delta
            avg_count+=1
            if avg_count == 5:
                avg_count = 0
                analyze_avg.append(0.0)
        del analyze_avg[-1]
        return analyze, np.array(analyze_avg)
    fid = open('pressure.log','r')
    P = fid.readlines()
    fid.close()
    fid = open('mylog.log','r')
    T = fid.readlines()
    fid.close()
    try:
        dname = os.getcwd().split('/')[-1].split('_')
        V = float(dname[dname.index('L')+1])**3
        N = int(dname[dname.index('nsphere')+1])
    except:
        fid = open('box_length.txt')
        fid.readline()
        S = fid.readline()
        fid.close()
        V = float(S.split()[1])**3
        N = int(S.split()[2])
    pressure, pressure_avg = log(P[1:],1)
    temperature, temperature_avg = log(T[1:],2)
    P =  sum(pressure_avg[xmin:xmax])/(pressure_avg.shape[0]-xmin-1)
    T =  sum(temperature_avg[xmin:xmax])/(temperature_avg.shape[0]-xmin-1)
    print "  Pressue:",P,"  Temperature:",T
    print "  N * T / V = ", N*T/ V
    pvt_avg[0].append(P)
    pvt_avg[1].append(N*T/V)
    plt.plot(pressure_avg[xmin:xmax], temperature_avg[xmin:xmax]*N/V, '%c'%colors[counter], label=title)
    counter+=1


if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
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
            print f
            os.chdir(f)
            title = f.split('_')[-1]
            pvt(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    index = int(config['-i']))
            os.chdir('../')
    else:
        dirs = []
        for i in sorted(os.listdir('.')):
            if os.path.isdir(i):
                dirs.append(i)
        for i in arg:
            f = dirs[int(i)]
            print f
            os.chdir(f)
            pvt(title = title,xmin = int(config['-x']), xmax = int(config['-r']),
                    index = int(config['-i']))
            os.chdir('../')
    plt.xlabel = 'Pressure'
    plt.ylabel = 'Temperature'
    plt.legend(loc=2)
    A = np.vstack([pvt_avg[0],np.ones(len(pvt_avg[0]))]).T
    r = np.linalg.lstsq(A,pvt_avg[1])[0]
    y = []
    for i in range(len(pvt_avg[0])):
        y.append(r[1]+r[0]*pvt_avg[0][i])
    print '#############################'
    print ' P = m * (N * T / V) + b, m = k'
    print '   m = ',r[0],'   b = ',r[1]
    plt.plot(pvt_avg[0],y,label='line of best fit')
    plt.show()


