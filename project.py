"Imports"

import numpy as np


"Données"
AR= 7
SW=10.76 #m²
bW=8.7 #m
c_root=1.9 #m
c_tip=0.57 #m
cVT_init = 0.08  #between 0.02-0.09
cHT_init = 0.75 #between 0.5-1
CW_bar =




"Calculations"
SVT_init = cVT_init*bW*SW/SVT
SHT_init = cHT_init*CW_bar*SW/LHT
dihedral_angle_init = np.arctan(np.sqrt(SVT_init/SHT_init))
