
import sys
import os
import numpy as np
import MD
############################################
## returns the arrays with the drift removed
############################################
def drift_remove_all(VW,V,W,L,step=5e4):
    from MD.analysis.drift_remove import eliminate_drift
    W = eliminate_drift(W,L)
    V = eliminate_drift(V,L)
    VW = eliminate_drift(VW,L)
    return VW,V,W
    #Find the msd of the system
if __name__ == '__main__':
#   #run_debug()
    #run_all()
    import MD
    print 'MD.L is '
    M=MD.ReadCord()
    Lx = M.box_length
    Ly = M.box_length_y
    Lz = M.box_length_z
    L = np.array([Lx, Ly, Lz])
    VW=M.cord_auto(['V','W'])
    drift_remove_all(VW,VW,VW,L)
