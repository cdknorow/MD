import os
import sys
sys.path.append('/home/cdknorow/Dropbox/Software/') 
import numpy as np
import math
import MD
import MD.plot.pyplot as pyplot
import MD.util as util


def scatter_plot_sort(M,Label,sort, delta=100):
    x = []
    y = []
    print Label
    count = 0
    yf = np.zeros((delta))
    new_label = []
    for j,i in enumerate(M):
        if Label[j] == sort:
            new_label.append(Label[j])
            if i[1][0] < 0.5:
                x.append(i[0])
                y.append(i[1]*5.)
            else:
                x.append(i[0])
                y.append(i[1])
            yf += y[count][0:delta]
            count+=1
    yf = yf/len(new_label)
    new_label.append('average')
    print len(x[0])
    x.append(x[0][0:delta])
    y.append(yf)
    pyplot.plot(x[0][0:delta],yf,'average',logx=True,limx=(0,100),limy=(0,1),save='average')
    Mark = ['x','x','o','o','v','v','1','1','8','8','D','D']
    Mark = ['x','x','o','o','v','1','1','8','8','D','D']
    pyplot.plot_multi(x,y,new_label,logx=True,limx=(0,100),limy=(0,1),loc=1,make_marker=True,M=Mark,
            showleg=True,save='kowasaki_%s'%sort)
    fid = open('intstat_avg_%s.dat'%sort,'w')
    for j in range(yf.shape[0]):
        fid.write(('%.2f %.4f\n'%(x[0][j],yf[j])))
    fid.close()

def scatter_plot(M,Label,delta=100):
    x = []
    y = []
    count = 0
    for i in M:
        if i[1][0] < 0.5:
            print 'Boosting', Label[count],count
            x.append(i[0])
            y.append(i[1]*5.)
        else:
            x.append(i[0])
            y.append(i[1])
        count+=1
    Mark = ['x','x','o','o','v','1','1','8','8','D','D']
    Mark = ['x','x','o','o','v','v','1','1','8','8','D','D','o','o']
    print len(x),len(Mark)
    pyplot.plot_multi(x,y,Label,logx=True,limx=(0,100),limy=(0,1),loc=1,make_marker=True,M=Mark,
            showleg=True,save='kowasaki')
    fid = open('intstat.dat','w')
    for i in Label:
        fid.write((' %s'%(i)))
    fid.write('\n')
    for i in range(delta):
        fid.write((' %.2f'%(i)))
        for j in range(len(y)):
            fid.write((' %.4f'%(y[j][i])))
        fid.write('\n')
    fid.close()
    #fid.write('time ')
    #for j in range(len(Label)):
    #    fid.write('%s '%Label[j].split())
    #fid.write('\n')
    #for i in range(len(x[0])):
    #    fid.write('%i '%x[0][i])
    #    for j in range(len(y)):
    #        fid.write('%.7f '%y[j][i])
    #    fid.write('\n')

if __name__ == '__main__':
    M = []
    Label = []
    delta = int(sys.argv[1])
    for f in sorted(os.listdir("./")):
        if os.path.isdir(f) and f[0]=='p':
            print f
            os.chdir(f)
            p = util.pickle_load('si_scat.pkl')
            if len(p[0]) > delta:
                M.append(p)
                Label.append('%s'%(f.split('_')[1:]))
            os.chdir('../')
        if os.path.isdir(f) and f[0]=='F':
            print f
            os.chdir(f)
            try:
                p = util.pickle_load('si_scat.pkl')
                print len(p[0])
                if len(p[0]) > delta:
                    M.append(p)
                    Label.append('%s'%(f.split('_')[9]))
            except:
                print 'no pickle'
            os.chdir('../')
    print Label
    if len(sys.argv)>2:
        delta = int(sys.argv[1])
        scatter_plot_sort(M,Label,sys.argv[2],delta)
    else:
        scatter_plot(M,Label,int(sys.argv[1]))
