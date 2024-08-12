# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

import numpy as np

ZL = np.zeros((5,20))

# iterate on height of the lattice of data qubits
for H_index,H in enumerate([4,5,6,7,8]):
    # iterate on physical error rate
    for p in range(20):
        try:
            data = open('./out/MonteCarlo_'+str(H)+'_'+str(p)+'.txt')
            lines = data.readlines()
            nb = len(lines) 
            ZL[H_index,p] = np.real(lines[nb-1])
        except:
            print('Could not read file data_'+str(H)+'_'+str(p)+'.txt')
            ZL[H_index,p] = np.nan
        data.close()
        
np.save("ZLoptimized_pheno",ZL)
