# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

# class to represent a column (= a timestep) in the circuit that contains gates
class Column:
    
    def __init__(self,nqubit):
        # number of qubits
        self.n = nqubit
        
        # list of gates         
        self.gates = []
        
        # number of gates in the list of gates
        self.nbgates = 0
        
    # add Gate to the column
    def AddGate(self,gate):
        # check if gate acts on a valid qubit
        for i in range(gate.size):
            if gate.qubits[i]>self.n-1:
                print("gate act on qubit out of bounds")
                
        self.gates.append(gate)
        self.nbgates+=1

    # propagate psi through the self column 
    def Propagate(self,psi,Mmt):
                
        for i in range(self.nbgates):
            
            gate = self.gates[i]
            
            # CNOT gate
            if gate.name == "CX":
                psi[gate.qubits[0]]=(psi[gate.qubits[0]]+psi[gate.qubits[1]])%2
                
            # X measurement
            if gate.name == "MX":
                # measurement error
                if gate.error == 1:
                    psi[gate.qubits[0]]=(psi[gate.qubits[0]]+1)%2
                # reset qubit
                if psi[gate.qubits[0]]==1:
                    psi[gate.qubits[0]]=0
                    Mmt.append(1)
                else:
                    Mmt.append(0)
                    
            # Preparation
            if gate.name == "P":
                psi[gate.qubits[0]] = 0
                                                        
    # trigger errors on gates              
    def Error(self,psi):
  
        for i in range(self.nbgates):
            
            gate = self.gates[i]

            # CNOT gate
            if gate.name == "CX" :
                
                if gate.error == 3 :
                    psi[gate.qubits[0]]=(psi[gate.qubits[0]]+1)%2
                    psi[gate.qubits[1]]=(psi[gate.qubits[1]]+1)%2
                    
                elif gate.error == 2:
                    psi[gate.qubits[1]]=(psi[gate.qubits[1]]+1)%2
                    
                elif gate.error == 1:
                    psi[gate.qubits[0]]=(psi[gate.qubits[0]]+1)%2
                    
            # Preparation
            if gate.name == "P":
                if gate.error == 1 :
                    psi[gate.qubits[0]]=(psi[gate.qubits[0]]+1)%2
                    
            # Idle gate
            if gate.name == "I":
                if gate.error == 1 :
                    psi[gate.qubits[0]]=(psi[gate.qubits[0]]+1)%2
                    

    def print(self):
        print(self)
        for i in range(self.nbgates):
            self.gates[i].print()

