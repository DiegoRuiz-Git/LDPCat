# Figure 5

The code in this folder estimates the logical error of the cellular automaton code for various numbers of iterations of the belief propagation decoding.

The folder `Circuit level parity check matrix` is designed to compute the parity check matrix that will be fed to the decoder, as well as the probability of each error mechanism and the final error on data qubits created by each error mechanism. It contains five Python files and `main.py` generates the matrices that are then stored in the folder `ParityCheck`.

The folder `Logical error` estimates the logical error rate of the BP+OSD decoder from https://github.com/quantumgizmos/bp_osd using the matrices from `Circuit level parity check matrix`. The Monte Carlo simulation is performed in `MonteCarlo.py`, which generates data in the `out` folder. `ReadData.py` then compiles the data into `ZL_bp_test.npy`, which are then plotted in the Jupyter noteebook `Plot_Figure5.ipynb`.