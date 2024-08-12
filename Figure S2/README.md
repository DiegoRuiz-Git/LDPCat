# Figure S2

Figure S2 evaluates the logical error rate of cellular automaton codes as a function of the figure of merit of cat qubits, $\kappa_1/\kappa_2$. 

Three different error models are examined: a phenomenological error model, a generic phase-flip circuit-level error model, and a cat qubit error model corresponding to Figure S2 (a), (b) and (c), respectively.



## Figure S2 a

The code evaluates the logical error rate under phenomenological noise for three different codes: the repetition code, the T-shape cellular automaton code, and the optimized cellular automaton code with different stabilizer shapes per row. This corresponds to the folders `repetition code`, `T-shape code` and `optimized code`.

In the folder `repetition code`, the Python file `MonteCarlo.py` performs a Monte Carlo simulation to evaluate the logical error of the repetition code. The BP+OSD decoder (https://github.com/quantumgizmos/bp_osd) is used, and the results are stored in the folder `out`. The Python file `ReadData.py` reads the data in the `out` folder and stores them in the numpy array `ZLrepetition_pheno.npy`.

In the folder `T-shape code`, the folder `Circuit level parity check matrix` is designed to compute the parity check matrix that will be fed to the decoder, as well as the probability of each error mechanism and the final error on data qubits created by each error mechanism. It contains five Python files, and `main.py` generates the matrices, which are then stored in the `ParityCheck` folder.

The folder `Logical error` estimates the logical error rate of the BP+OSD decoder from https://github.com/quantumgizmos/bp_osd fed with the matrices from `Circuit level parity check matrix`. The Monte Carlo simulation is performed in `MonteCarlo.py`, which generates data in the `out` folder. `ReadData.py` then compiles the data into `ZLTshape_pheno.npy`.

The folder `optimized code` is organized in the same way as the `T-shape code` but computes the logical error of the optimized cellular automaton codes with different stabilizer shapes. \
Note: The `ParityCheck` folder has been zipped to save space for the optimized codes in Figures S2 a,b and c.

The Jupyter notebook `Plot_FigureS2a.ipynb` compiles the result of the three different codes and plots them.

## Figure S2 b

The folder for Figure S2 b is organized in the exact same way as Figure S2 a, except a circuit-level error model is used instead of a phenomenological error model.

## Figure S2 c

The folder of Figure S2 c is organized in the same way as Figure S2 a and b, except a noise model specific to the cat qubit architecture is used, and the simulation is only performed for the repetition and optimized codes.