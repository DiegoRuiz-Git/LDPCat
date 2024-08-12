# Figure 1

The code in this folder is designed to compute the parameters of local codes whose stabilizers have a range covering a patch of 3 $\times$ 3 data qubits. The local codes being investigated also possess a single stabilizer shape, which is then translated horizontally and vertically across the lattice to generate the other stabilizers.

## Figure 1a

In Figure 1a, we test all possible stabilizer shapes (each one corresponding to a code) with the above criteria on lattices of data qubits with sizes up to $17 \times 17$. For each code, we computed the number of logical qubits $k$, and the distance $d$.

When the number of logical qubits is smaller than 34, we compute the Hamming weight of the $2^k$ codewords to obtain the distance. If it is larger, the previous method does not work, and we resort to a SAT solver.

The folder `Figure 1a` contains a C++ program, `Calculate_code_parameters.cpp` which computes the code parameters and outputs the result in the `result` folder. The Jupyter notebook `Plot_Figure1a.ipynb` then analyses and plots the results. \
Note: The z3 library is required to compile and run the code. You can find it at https://github.com/Z3Prover/z3.

## Figure 1b

The code `Calculate_code_parameters.cpp` for Figure 1b is very similar, the only difference being that we only consider codes with a number of logical qubits $k = 22$ and that we look for codes on lattices of data qubits of sizes up to $34 \times 34$. As with Figure 1a, the results are saved in the `result` folder and analyzed and plotted with the jupyter notebook `Plot_Figure1b.ipynb`

The folder also contains a file `BestCodes.csv`, containing the data of the best found and possible classical codes from http://www.codetables.de/BKLC/Tables.php?q=2&n0=1&n1=256&k0=1&k1=256 that is also used for the plot.