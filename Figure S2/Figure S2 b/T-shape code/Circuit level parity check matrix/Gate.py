# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

# class to represent a gate in the circuit
class Gate:
    
    def __init__(self,name,size,qubits):
        # name of the gate
        self.name = name
        # number of qubits of the gate
        self.size = size
        # list of qubits
        self.qubits = qubits
        # should trigger error
        self.error = 0
        
    def print(self):
        print(self.name+" acts on "+str(self.size)+" qubit : ")
        print(self.qubits)
        