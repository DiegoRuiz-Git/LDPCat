# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

from bposd.css import css_code
from ldpc import bposd_decoder
import numpy as np
import os
import time
import sys
import logging
import multiprocessing

logging.basicConfig(
    level=logging.INFO,
)

# number of cores for multiprocessing
num_cores = multiprocessing.cpu_count()
print("Number of CPU cores:", num_cores)

############################### PREPARE DECODER ############################### 

# distance of the repetition code
d = int(sys.argv[1])

# create the circuit level parity check matrix of the repetition code
def ParityMatrix(d):
    
    # data errors
    H = np.zeros((d-1,d),dtype = int)
    for i in range(d-1):
        line = np.zeros((d))
        line[0] = 1
        line[1] = 1
        H[i] = np.roll(line,i)
        
    A = np.kron(np.eye(d+1,dtype = int),H)

    # measurement errors
    B = np.zeros(((d-1)*(d+1),(d-1)*d),dtype = int)
    for i in range(len(B[0])):
        column = np.zeros(((d-1)*(d+1),1))
        column[0] = 1
        column[d-1] = 1
        B[:,i] = np.roll(column[:,0],shift=i)
        
    return np.hstack((A,B))

# physical error rate on data qubits
p = np.logspace(-3,np.log10(0.3),20)[int(sys.argv[2])]
# physical error rate on measurements
q = p

# Probability matrix
def ProbaMatrix(d):
    
    # data qubit errors
    A = np.full((d*(d+1)),p)
                
    # measurement errors
    B = np.full(((d-1)*d),q)
                    
    return np.hstack((A,B))

H = ParityMatrix(d)
ProbaMatrix = ProbaMatrix(d)
# create the repetition code
rep_code = css_code(hx = H,hz = np.zeros((1,len(H[0])),dtype = int)) 
        
# initialize the decoder
bpd=bposd_decoder(
    rep_code.hx, #the parity check matrix
    max_iter = 10000, #the maximum number of iterations for BP
    bp_method = "ms", #The BP method. using "ms": min-sum updates
    ms_scaling_factor = 0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    channel_probs = ProbaMatrix, #assign error_rate to each qubit
    osd_method="osd_cs", #the OSD method. Choose from: "osd_e" (Panteleev and Kalechev), "osd_cs" (Roffe), "osd0" (Order 0 with e_[T] = 0)
    osd_order = min(len(H[0]) - len(H),60) #the lambda in Roffe paper <= nb_edge - nb_syndrom
)

############################### MONTECARLO ############################### 

def MonteCarlo():
  
    error_count = 0
    run_count = 0

    # generate seed
    seed = (int(time.time() * 1000) + os.getpid()) % (2**32)
    np.random.seed(seed)
    print(f"Seed for process {os.getpid()}: {seed}")

    # run for 2 errors for each core
    while error_count< 2:

        # Generate error
        random_array = np.random.rand(len(H[0]))
        error = (random_array < ProbaMatrix).astype(int)
    
        # compute syndromes
        syndrome = rep_code.hx@error %2

        # decode
        bpd.decode(syndrome)

        # error after correction
        residual_error=(bpd.osdw_decoding+error) %2

        # final error on data qubits
        final_error = np.zeros((d),dtype = int)
        for j in range(len(final_error)):
            for k in range(d+1):
                final_error[j] += residual_error[k*d+j]

        final_error = final_error%2 

        run_count+=1

        if not np.all(final_error == 0):
            error_count += 1
            logging.info("error_count "+str(error_count))

    return run_count

############################### RUN WITH MULTIPROCESSING ############################### 

with multiprocessing.Pool() as pool:
    results = [pool.apply_async(MonteCarlo) for _ in range(64)]
    # Retrieve the results
    results = [result.get() for result in results]

run_count = 0
for i in range(len(results)):
    run_count += results[i]

print(128/run_count)