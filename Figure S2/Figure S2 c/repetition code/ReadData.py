# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

# put result from Monte Carlo in numpy array
import numpy as np

ZL = np.zeros((5,20))

# iterate on distance
for d_index,d in enumerate([5,9,12,16,22]):
    # iterate on physical error rate
    for p in range(20):
        data = open('./out/data_'+str(d)+'_'+str(p)+'.txt')
        lines = data.readlines()
        nb = len(lines) 
        try:
            ZL[d_index,p] = np.real(lines[nb-1])
        except:
            print('Could not read file data_'+str(d)+'_'+str(p)+'.txt')
            ZL[d_index,p] = np.nan
        data.close()
        
np.save("ZLrepetition_cat",ZL)
