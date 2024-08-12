# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

# put result from Monte Carlo in numpy array
import numpy as np

# parameters in order:
# belief propagation method
# nb iterations bp
# scaling factor
# physical error
ZL = np.zeros((4,5,4,13))

for m in range(4):
    for r_index,r in enumerate([4,1032,5247,10000,11365]):
        for s in range(4):
            for p in range(13):
                data = open('./out/MonteCarlo_'+str(m)+'_'+str(r)+'_'+str(s)+'_'+str(p)+'.txt')
                lines = data.readlines()
                nb = len(lines) 
                ZL[m,r_index,s,p] = np.real(lines[nb-1])
                data.close()
        
np.save("ZL_bp_test",ZL)
