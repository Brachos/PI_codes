"Imports"

import numpy as np
from scipy.integrate import quad

"Données"
AR= 7
SW=10.76 #m²
bW=8.7/2 #m
c_root=1.9 #m
c_tip=0.57 #m
cVT_init = 0.08  #between 0.02-0.09
cHT_init = 0.75 #between 0.5-1
LVT = 1.5
LHT = 3

"Calculations"
def integrand(x, c_root, c_tip, bW):
    return (c_root + (c_tip-c_root)/bW*x)**2
CW_bar = 2/SW * quad(integrand,0,bW,args=(c_tip,c_root,bW))
SVT_init = cVT_init*bW*SW/LVT
SHT_init = cHT_init*CW_bar*SW/LHT
dihedral_angle_init = np.arctan(np.sqrt(SVT_init/SHT_init))
print("L'angle initial dihedre est de %.3f degrés" %(dihedral_angle_initcd))
