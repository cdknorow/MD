
import sys
import os
import pickle
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#set the fonts
matplotlib.rc('font', **{'size':'20'})
plt.close()
plt.figure(figsize=(20,20), dpi=100)
colors = ['b','r','g','b','k','y','c','m']
counter = 0

fig1 = plt.figure(figsize=(20,8), dpi=100)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

############################################
def function(index=1):
    global counter
    global colors
    log = open('mylog.log','r')
    title = log.readline().split()[index]
    logger = [[],[]]
    time = []
    for line in log.readlines():
        logger[0].append(float(line.split()[2]))
        logger[1].append(float(line.split()[3]))
        time.append(line.split()[0])
    log.close()
    Energy = np.array(logger[0])+np.array(logger[1])
    DE = np.array(logger[0]) - np.array(logger[1])
    ax1.plot(time,Energy,'b', lw=.5, label = 'Total Energy')
    ax2.plot(time,DE,lw=.5, label = 'KE - PE')
    ax2.set_xlabel('time')
    ax1.set_ylabel('Totale Energy')
    ax2.set_ylabel('KE-PE')
    plt.show()

if __name__ == "__main__":
    import getopt
    # -x max xrange for fitting function
    # -f filename
    # -r max xrange for plotting
    # arg 1 2 3 4,etc directory number in current dir
    config = {'-x':'-1','-r':'-1','-i':0}
    options, arg =  getopt.getopt(sys.argv[1:],'x:r:D:i:')
    config.update( dict(options))
    function(index = int(config['-i']))
    plt.show()

