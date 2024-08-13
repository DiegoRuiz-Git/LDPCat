# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

from bposd.css import css_code
from ldpc import bposd_decoder
import numpy as np
import math
import os
import sys
import logging
import multiprocessing
import time

logging.basicConfig(
    level=logging.INFO,
)

# number of cores for multiprocessing
num_cores = multiprocessing.cpu_count()
logging.info("Number of CPU cores:"+ str(num_cores))

############################### PREPARE DECODER ############################### 

# store the triggered syndromes by every error mechanism
H = np.load("../Circuit level parity check matrix/ParityCheck/CircuitLevelMatrix" + sys.argv[1] + '_' + sys.argv[2] + ".npy")
# store the final error on data qubits from every error mechanism
ErreurFinale = np.load("../Circuit level parity check matrix/ParityCheck/ErrorMatrix" + sys.argv[1] + '_' + sys.argv[2] + ".npy")
# store the probability of every error mechanism
ProbaMatrix = np.load("../Circuit level parity check matrix/ParityCheck/ProbaMatrix" + sys.argv[1] + '_' + sys.argv[2] + ".npy")

logging.info("H dim :" + str(H.shape))

# create the cellular automaton code
optimized_code = css_code(hx = H,hz = np.zeros((1,len(H[0])),dtype = int)) 

# initialize the decoder
bpd=bposd_decoder(
    optimized_code.hx, #the parity check matrix
    max_iter = 10000, #the maximum number of iterations for BP
    bp_method = "ms", #The BP method. using "ms": min-sum updates
    ms_scaling_factor = 0.625, #min sum scaling factor. If set to zero the variable scaling factor method is used
    channel_probs = ProbaMatrix, #assign error_rate to each qubit.
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
    # more errors were sampled for high physical error rate
    while error_count < 2:

        # Generate error
        random_array = np.random.rand(len(H[0]))
        error = (random_array < ProbaMatrix).astype(int)

        # Real final error
        indices = np.where(error == 1)[0]

        RealFinalError = np.zeros((len(ErreurFinale)),dtype = int)
        for i in indices:
            RealFinalError = (RealFinalError + ErreurFinale[:,i].T)%2

        # compute syndromes
        syndrome = optimized_code.hx@error.T %2

        # decode
        bpd.decode(syndrome)

        # compute correction
        indices = np.where(bpd.osdw_decoding == 1)[0]

        correction = np.zeros((len(ErreurFinale)),dtype = int)
        for i in indices:
            correction = (correction + ErreurFinale[:,i].T)%2

        # compare correction and the real final error
        final_error = (correction + RealFinalError)%2

        run_count +=1

        if not np.all(final_error == 0):
            error_count+=1
            logging.info("error_count : "+str(error_count))

    return error_count,run_count

############################### RUN WITH MULTIPROCESSING ###############################

if __name__ == '__main__':

    with multiprocessing.Pool() as pool:
        results = [pool.apply_async(MonteCarlo) for _ in range(50)]
        # Retrieve the results
        results = [result.get() for result in results]

    error_count = 0
    run_count = 0
    for i in range(len(results)):
        error_count += results[i][0]
        run_count += results[i][1]

    print(error_count/run_count)
