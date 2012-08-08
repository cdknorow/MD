import numpy as np
import MD
import MD.analysis.bond_order as bond_order 
from MD.util import util as util

reload(bond_order)
VW = util.pickle_load('bcc_unit.pkl')
L=4
#bond_order.structure_correlation(VW,L,l=6,crystal=8,rcut=0.88,verbose=False)
bond_order.solid_particles(VW,L,l=6,crystal=8,rcut=0.88,verbose=False)

#bond_order.steinhardt_order(VW,L,l=2,rcut=1.01,verbose=False)
#bond_order.steinhardt_order(VW,L,l=4,rcut=1.01,verbose=False)
#bond_order.steinhardt_order(VW,L,l=6,rcut=1.01,verbose=False)
#bond_order.steinhardt_order(VW,L,l=8,rcut=1.01,verbose=False)
#bond_order.steinhardt_order(VW,L,l=10,rcut=1.01,verbose=False)

