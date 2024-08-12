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

############################### PREPARE DECODER ############################### 

# distance of the repetition code
d = int(sys.argv[1])

# Cat qubit parameters
nbar = 11 # number of photons
k1 = np.logspace(-5,np.log10(3e-3),20)[int(sys.argv[2])] # single photon loss rate
k2 = 1 # two-photon dissipation rate

# create the circuit level parity check matrix of the repetition code
def ParityMatrix(d):
    
    # space errors
    H = np.zeros((d-1,d),dtype = int)
    for i in range(d-1):
        line = np.zeros((d))
        line[0] = 1
        line[1] = 1
        H[i] = np.roll(line,i)
        
    A = np.kron(np.eye(d+1,dtype = int),H)

    # time errors
    B = np.zeros(((d-1)*(d+1),(d-1)*d),dtype = int)
    for i in range(len(B[0])):
        column = np.zeros(((d-1)*(d+1),1))
        column[0] = 1
        column[d-1] = 1
        B[:,i] = np.roll(column[:,0],shift=i)

    # space-time errors
    C = np.zeros(((d-1)*(d+1),d*(d-2)),dtype = int)
    for i in range(d):
        for j in range(d-2):
            column = np.zeros(((d-1)*(d+1),1))
            column[(d-1)*i] = 1
            column[1 + (d-1)*(i+1)] = 1
            C[:,i*(d-2) + j] = np.roll(column[:,0],shift=j)
        
    return np.hstack((A,B,C))

# output the probability that an error occurred by combining a list of probabilities
def combine(list_p):
    if len(list_p) == 1:
        return list_p[0]
    if len(list_p) == 2:
        return list_p[0]*(1-list_p[1]) + (1-list_p[0])*list_p[1]
    else:
        p=combine(list_p[1:])
    return list_p[0]*(1-p)+(1-list_p[0])*p

# physical error rate of error mechanisms
p_prep = nbar*k1/k2 # preparation error
p_meas = nbar*k1/k2 # measurement error
p_CX_Z1 = nbar*k1/k2 + np.pi**2/64/nbar # CNOT ancilla error
p_CX_Z2 = 1/2*nbar*k1/k2 # CNOT target error
p_CX_Z1Z2 = 1/2*nbar*k1/k2 # CNOT correlated error
p_refresh = nbar*k1/k2 # idle error during wait time to remove leakage on ancilla qubit after CNOT
p_prep_w = nbar*k1/k2 # Idle error during preparation 
p_meas_w = nbar*k1/k2 # Idle error during measurement 
p_CX_w = nbar*k1/k2 # Idle error during CNOT 
p_refresh_w = nbar*k1/k2 # idle error during wait time to remove leakage on data qubit after CNOT

# probabilities associated to edges of the error graph
# horizontal edge (measurement errors)
p_hor = combine([p_prep,p_CX_Z1,p_refresh,p_CX_Z1,p_meas])
# vertical edge (data errors)
p_ver = combine([p_CX_Z2,p_meas_w,p_prep_w,p_CX_Z1Z2])
# input vertical edge
p_in = combine([p_prep_w,p_CX_Z1Z2])
# input bottom vertical edge
p_in_bottom = combine([p_prep_w,p_CX_w,p_refresh_w,p_CX_Z1Z2])
# vertical edge on the first row
p_top = combine([p_CX_Z2,p_refresh_w,p_CX_w,p_meas_w,p_prep_w,p_CX_Z1Z2])
# vertical edge on the last row
p_bottom = combine([p_CX_Z2,p_meas_w,p_prep_w,p_CX_w,p_refresh_w,p_CX_Z1Z2])
# output vertical edge
p_out = combine([p_CX_Z2,p_meas_w])
# top output vertical edge
p_out_top = combine([p_CX_Z2,p_refresh_w,p_CX_w,p_meas_w])
# diagonal edge
p_diag = combine([p_CX_Z2,p_refresh_w,p_CX_Z1Z2])

# Probability matrix
def ProbaMatrix(d):
    
    # input vertical edges
    A = np.zeros((d))
    for i in range(d-1):
        A[i] = p_in
    A[d-1] = p_in_bottom
    
    # vertical edges
    B = np.zeros((d))
    B[0] = p_top
    for i in range(1,d-1):
        B[i] = p_ver
    B[d-1] = p_bottom
    
    # output vertical edges
    C = np.zeros((d))
    C[0] = p_out_top
    for i in range(1,d):
        C[i] = p_out
        
    # stack vertical edges
    D = A
    for _ in range(d-1):
        D = np.hstack((D,B))
    D = np.hstack((D,C))
            
    # horizontal edges
    E = np.zeros(((d-1)*d))
    for i in range(len(E)):
        E[i] = p_hor
        
    # diagonal edges
    F = np.zeros(((d-2)*d))
    for i in range(len(F)):
        F[i] = p_diag
            
    return np.hstack((D,E,F))

H = ParityMatrix(d)
ProbaMatrix = ProbaMatrix(d)
# create the repetition code
rep_code = css_code(hx = H,hz = np.zeros((1,len(H[0])),dtype = int)) 
        
bpd = bposd_decoder(
    rep_code.hx, #the parity check matrix
    max_iter = 10000, #the maximum number of iterations for BP
    bp_method = "ms", #The BP method. using "ms": min-sum updates
    ms_scaling_factor = 0.625, #min sum scaling factor. If set to zero the variable scaling factor method is used
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
    while error_count < 2:

        # Generate error
        random_array = np.random.rand(len(H[0]))
        error = (random_array < ProbaMatrix ).astype(int)
    
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

        for j in range(len(final_error)-2):
            for k in range(d):
                final_error[j+1] += residual_error[d*(d+1) + (d-1)*d + (d-2)*k + j]

        final_error = final_error%2 

        run_count+=1

        if not np.all(final_error == 0):
            error_count+=1
            logging.info("error_count : "+str(error_count))

    return error_count,run_count

############################### RUN WITH MULTIPROCESSING ############################### 

with multiprocessing.Pool() as pool:
    results = [pool.apply_async(MonteCarlo) for _ in range(64)]
    # Retrieve the results
    results = [result.get() for result in results]

error_count = 0
run_count = 0
for i in range(len(results)):
    error_count += results[i][0]
    run_count += results[i][1]

print(error_count/run_count)