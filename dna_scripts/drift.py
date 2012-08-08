
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
    print MD.L
    MD.V = MD.M.cord_auto(['V'])
    MD.W = MD.M.cord_auto(['W'])
    MD.VW = MD.M.cord_auto(['V','W'])
    drift_remove_all(MD.VW,MD.V,MD.W,MD.L)
