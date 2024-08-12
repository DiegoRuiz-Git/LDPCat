# Table II

The code in this folder is designed to find the shapes of stabilizers that maximize the distance of cellular automaton codes. This folder contains two C++ programs `BestShapes_Periodic_SAT.cpp` and `BestShapes_NonPeriodic_SAT.cpp`.

The `BestShapes_Periodic_SAT.cpp` implements a cellular automaton code with periodic boundary conditions on the lateral sides, while `BestShapes_NonPeriodic_SAT.cpp` implements a cellular automaton code without periodic boundary conditions as shown in Figure S1.

Both programs use the z3 SAT solver (https://github.com/Z3Prover/z3) to find the maximal attainable distance and the corresponding stabilizer shapes. For a given lattice of data qubits of height $H$ and width $L$, the SAT solver checks whether there is a set of stabilizer shapes that achieve a certain distance $d$. If the SAT solver returns that this problem is impossible, it means this distance is out of reach of the cellular automaton codes we are exploring.
 
The non-periodic case is the one that will ultimately be implemented experimentally and is thus the primary focus of our interest. However, increasing the number of logical qubits can decrease (but never increase) the distance due to interference effects. This program can thus provide an upper bound on the distance of cellular automaton codes for a large number of logical qubits.
 
On the other hand, in the periodic case, the distance can increase (but never decrease) as the number of logical qubits increases. Indeed, due to the periodic boundary conditions, some codewords can interfere destructively with themselves by 'looping' around the lattice if the number of logical qubits is too small. Thus the distance given by this program represents a lower bound of the distance of the non-periodic case for any number of logical qubits.

Thus, by varying $L$ for various $H$ for both programs, we are able to derive an upper bound and a lower bound for the distance for each $H$. We observed empirically that this lower bound and upper bound are always equal if we go to large enough $L$.

