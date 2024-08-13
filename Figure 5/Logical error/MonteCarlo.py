# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

from bposd.css import css_code
from ldpc import bposd_decoder
import numpy as np
import os
import sys
import math
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
H = np.load("../Circuit level parity check matrix/ParityCheck/CircuitLevelMatrix" + "3" + '_' + sys.argv[4] + ".npy")
# store the final error on data qubits from every error mechanism
FinalError = np.load("../Circuit level parity check matrix/ParityCheck/ErrorMatrix" + "3" + '_' + sys.argv[4] + ".npy")
# store the probability of every error mechanism
ProbaMatrix = np.load("../Circuit level parity check matrix/ParityCheck/ProbaMatrix" + "3" + '_' + sys.argv[4] + ".npy")

logging.info("H dim :" + str(H.shape))

# create the cellular automaton code
tetris_code = css_code(hx = H,hz = np.zeros((1,len(H[0])),dtype = int)) 

# number of iterations of BP
bp_rep = int(sys.argv[2])
logging.info("bp rep : " + str(bp_rep))

# scaling factor for BP
scaling_factor = 0
if int(sys.argv[3])==0:
    scaling_factor = 0
if int(sys.argv[3])==1:
    scaling_factor = 0.1
if int(sys.argv[3])==2:
    scaling_factor = 0.625
if int(sys.argv[3])==3:
    scaling_factor = 1

logging.info("scaling_factor : " + str(scaling_factor))

# BP method
bp_method = 0
if int(sys.argv[1])==0:
    bp_method = "ps"
if int(sys.argv[1])==1:
    bp_method = "ms"
if int(sys.argv[1])==2:
    bp_method = "msl"
if int(sys.argv[1])==3:
    bp_method = "psl"

logging.info("bp_method : " + str(bp_method))

# initialize the decoder
bpd=bposd_decoder(
    tetris_code.hx, #the parity check matrix
    max_iter = bp_rep, #the maximum number of iterations for BP
    bp_method = bp_method, #The BP method
    ms_scaling_factor = scaling_factor, #min sum scaling factor. If set to zero the variable scaling factor method is used
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

    # run for 100/num_cores errors for each core
    while error_count < math.ceil(100/num_cores):

        # Generate error
        random_array = np.random.rand(len(H[0]))
        error = (random_array < ProbaMatrix).astype(int)

        # Real final error
        indices = np.where(error == 1)[0]

        RealFinalError = np.zeros((len(FinalError)),dtype = int)
        for i in indices:
            RealFinalError = (RealFinalError + FinalError[:,i].T)%2

        # compute syndromes
        syndrome = tetris_code.hx@error.T %2

        # decode
        bpd.decode(syndrome)

        # compute correction
        indices = np.where(bpd.osdw_decoding == 1)[0]

        correction = np.zeros((len(FinalError)),dtype = int)
        for i in indices:
            correction = (correction + FinalError[:,i].T)%2

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