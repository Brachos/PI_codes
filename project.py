"Imports"

import numpy as np


"Données"
#Gamma #Dihedral of the V-tail
#Cl_alpha_N #lift curve slope of the planform at a zero dihedral
#tau #control effectiveness parameter
#K #the ratio of sum of lifts obtained by equal and
#opposite changes in angle-of-attack of two semispans
#of tail to lift obtained by an equal change in
#angle-of-attack for complete tail

"Lift curve slope of the V-tail"
Cl_alpha_HT = Cl_alpha_N*(np.cos(Gamma))**2
"Elevator authority of the V-tail"
Cl_delta_e = Cl_alpha_N*tau*np.cos(Gamma)
"Side force slope of the V-tail"
CY_beta_r = -K*Cl_alpha_N*(np.sin(Gamma))**2
"Rudder authority of the V-tail"
CY_delta_r = K*Cl_alpha_N*np.sin(Gamma)
