# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

from FractalCode import FractalCode
import numpy as np
import sys

H = int(sys.argv[1]) # height of the lattice of data qubits
L = 10 # width of the lattice of data qubits
distance = [7,12,15][H-3] # distance of the code

# physical error rate
p = np.logspace(-3,np.log10(0.3),20)[int(sys.argv[2])]

fractalcode = FractalCode(H,L)

fractalcode.CreateCheckMatrix()

# number of rounds = distance
fractalcode.CreateCircuit(distance)

print("circuit depth ",fractalcode.circuit.nbcol)

# calculate number of error mechanisms
nb_pauli_errors = 0
for i in range(fractalcode.circuit.nbcol):
    for j in range(fractalcode.circuit.columns[i].nbgates):
        if fractalcode.circuit.columns[i].gates[j].name == "CX":
            nb_pauli_errors += 3
        else:
            nb_pauli_errors += 1

# size of the circuit level parity check matrix
print("number of space time syndromes ", (H-1)*L*(distance+1) )
print("number of error mechanisms ",nb_pauli_errors)

############## Compute the circuit level parity check matrix ##############

# store the triggered syndromes by every error mechanism
CircuitLevelMatrix = np.zeros(((H-1)*L*(distance+1),nb_pauli_errors))
# store the final error on data qubits from every error mechanism
ErrorMatrix = np.zeros((H*L,nb_pauli_errors))
# store the probability of every error mechanism
ProbaMatrix = np.zeros((nb_pauli_errors))

err_counter = 0

## loop through columns
for i in range(fractalcode.circuit.nbcol):
    
    print(i," / ",fractalcode.circuit.nbcol)
    column = fractalcode.circuit.columns[i]
    
    ## loop through gates
    for j in range(column.nbgates):
        
        gate = fractalcode.circuit.columns[i].gates[j]
        
        ## loop through possible error on the gate
        possible_error = int(2*gate.size - 1)
        for k in range(possible_error):
                
                # trigger error on gate
                gate.error = k+1
                fractalcode.circuit.Error(bprint = True)
                
                # fetch syndrome measurement
                CircuitLevelMatrix[:(H-1)*L*distance,err_counter] = fractalcode.circuit.Mmt
                
                # fetch final error
                final_error = fractalcode.circuit.psi[:H*L]
                ErrorMatrix[:,err_counter] = final_error
                
                # compute last round of perfect syndrome measurement
                for l in range(len(fractalcode.CheckMatrix)):
                    syndrome = 0
                    for m in range(len(fractalcode.CheckMatrix[l])):
                        syndrome += final_error[fractalcode.CheckMatrix[l][m]]
                    CircuitLevelMatrix[(H-1)*L*distance + l,err_counter] = syndrome%2
                    
		        ## Compute error probability
                if gate.name == "CX":
                    if gate.error == 1:
                        ProbaMatrix[err_counter] = p/3
                    if gate.error == 2:
                        ProbaMatrix[err_counter] = p/3
                    if gate.error == 3:
                        ProbaMatrix[err_counter] = p/3
                if gate.name == "I":
                    ProbaMatrix[err_counter] = p
                if gate.name == "MX":
                    ProbaMatrix[err_counter] = p
                if gate.name == "P":
                    ProbaMatrix[err_counter] = p
                        
                # reset circuit
                fractalcode.circuit.Reset()
                gate.error = 0
                
                err_counter+=1
                
############## Convert measurement to detector  ##############

newCircuitLevelMatrix = np.zeros(((H-1)*L*(distance+1),nb_pauli_errors))

for l in range(len(CircuitLevelMatrix[0])):
    for n in range(L*(H-1)):
        newCircuitLevelMatrix[n,l] = CircuitLevelMatrix[n,l] 

# compute syndrome difference in time
for l in range(len(CircuitLevelMatrix[0])):
    for m in range(1,distance+1):
        for n in range(L*(H-1)):
            newCircuitLevelMatrix[m*L*(H-1) + n,l] = (CircuitLevelMatrix[m*L*(H-1) + n,l] 
                                                + CircuitLevelMatrix[(m-1)*L*(H-1) + n,l])%2
            
CircuitLevelMatrix = newCircuitLevelMatrix

############## Factorize the columns of the parity check matrix  ##############

StackMatrix = np.vstack((CircuitLevelMatrix, ErrorMatrix))

unique_columns = {}
new_StackMatrix = []
new_ProbaMatrix = []
new_StackMatrix_tuples = []

# Iterate through columns of StackMatrix
for i in range(StackMatrix.shape[1]):
    
    col = tuple(StackMatrix[:, i])
    
    # if equivalent error mechanism, sum probabilities
    if col in unique_columns:
        # Get index of the existing unique column
        index = new_StackMatrix_tuples.index(col)
        
        # Calculate the updated probability
        p1 = ProbaMatrix[i]
        p2 = new_ProbaMatrix[index]
        new_value = p1 * (1 - p2) + p2 * (1 - p1)
        
        # Update the probability in new_ProbaMatrix
        new_ProbaMatrix[index] = new_value
        
    # otherwise, add new column
    else:
        unique_columns[col] = 1
        new_StackMatrix.append(StackMatrix[:, i])
        new_StackMatrix_tuples.append(col)
        new_ProbaMatrix.append(ProbaMatrix[i])

# create final matrices
ProbaMatrix_final = np.array(new_ProbaMatrix)
StackMatrix_final = np.column_stack(new_StackMatrix)
CircuitLevelMatrix_final = StackMatrix_final[:(H-1)*L*(distance+1),:]
ErrorMatrix_final = StackMatrix_final[(H-1)*L*(distance+1):,:]

print("number of columns after factorization ",len(CircuitLevelMatrix_final[0]))

np.save("ParityCheck/CircuitLevelMatrix"+str(H)+"_"+sys.argv[2]+".npy",CircuitLevelMatrix_final)
np.save("ParityCheck/ErrorMatrix"+str(H)+"_"+sys.argv[2]+".npy",ErrorMatrix_final)
np.save("ParityCheck/ProbaMatrix"+str(H)+"_"+sys.argv[2]+".npy",ProbaMatrix_final)







