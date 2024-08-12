# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

import Circuit
import numpy as np
import Column
from Gate import Gate

# class that represents cellular automaton/fractal codes
class FractalCode:
    
    def __init__(self,H,L):
        # height of the code
        self.H = H
        # length of the code
        self.L = L

    # stabilizer for each row depending on H
    def Stabilizer(self):

        Stab = np.zeros((self.H-2,2,3),dtype = int)

        if self.H == 4:
            Stab[0] = np.array([[0,1,1],
                                [0,0,1]])
            Stab[1] = np.array([[1,0,0],
                                [0,1,1]])
        if self.H == 5:
            Stab[0] = np.array([[1,0,1],
                                [1,0,0]])
            Stab[1] = np.array([[1,0,1],
                                [1,0,0]])
            Stab[2] = np.array([[0,0,1],
                                [1,0,1]])

        if self.H == 6:
            Stab[0] = np.array([[1,0,1],
                                [0,0,1]])
            Stab[1] = np.array([[1,0,1],
                                [0,0,1]])
            Stab[2] = np.array([[1,0,0],
                                [1,0,1]])
            Stab[3] = np.array([[0,0,0],
                                [1,1,1]])

        if self.H == 7:
            Stab[0] = np.array([[1,0,1],
                                [0,1,0]])
            Stab[1] = np.array([[1,0,1],
                                [0,0,1]])
            Stab[2] = np.array([[0,0,1],
                                [1,0,1]])
            Stab[3] = np.array([[0,1,0],
                                [1,0,1]])
            Stab[4] = np.array([[0,0,0],
                                [1,1,1]])

        if self.H == 8:
            Stab[0] = np.array([[1,1,1],
                                [0,0,0]])
            Stab[1] = np.array([[1,1,1],
                                [0,0,0]])
            Stab[2] = np.array([[1,0,1],
                                [1,0,0]])
            Stab[3] = np.array([[1,0,1],
                                [1,0,0]])
            Stab[4] = np.array([[0,0,1],
                                [1,0,1]])
            Stab[5] = np.array([[1,0,0],
                                [1,0,1]])

        self.Stab = Stab

    # Parity check matrix of the code
    def CreateCheckMatrix(self):

        H = self.H
        L = self.L
        self.Stabilizer()
        Stab = self.Stab
       
        self.CheckMatrix = np.zeros(((H-2)*L,4),dtype = int )
               
        for i in range(H-2):
            for j in range(L):
                non_zero = np.where(Stab[i])
                self.CheckMatrix[i*L+j] = np.array([(i+2)*L+j,
                                                    (i + non_zero[0][0])*L + (j-1+non_zero[1][0])%L,
                                                    (i + non_zero[0][1])*L + (j-1+non_zero[1][1])%L,
                                                    (i + non_zero[0][2])*L + (j-1+non_zero[1][2])%L])

    # create the circuit of the cellular automaton code from the parity check matrix                
    def CreateCircuit(self,d):
        
        n = self.H*self.L # nb of data
        a = (self.H-2)*self.L # nb of ancilla
        
        self.circuit = Circuit.Circuit(n+a)
        
        for _ in range(d):
            
            # prepare ancilla
            column1 = Column.Column(n+a)
            for i in range(n,n + a):
                column1.AddGate(Gate("P",1,[i]))
                

            #idle on data
            for i in range(n):
                column1.AddGate(Gate("I",1,[i]))
            self.circuit.AddColumn(column1)
            

            # upper cnot
            column2 = Column.Column(n+a)
            for i in range(a):
                column2.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][0]]))

            # idle on qubits not touched
            all_data = set(range(n))
            for i in range(len(self.CheckMatrix)):
                all_data.discard(self.CheckMatrix[i][0])
            idle_data = list(all_data)
            for i in idle_data:
                column2.AddGate(Gate("I",1,[i]))
                
            self.circuit.AddColumn(column2)


            # right cnot
            column3 = Column.Column(n+a)
            for i in range(a):
                column3.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][1]]))

            # idle on qubits not touched
            all_data = set(range(n))
            for i in range(len(self.CheckMatrix)):
                all_data.discard(self.CheckMatrix[i][1])
            idle_data = list(all_data)
            for i in idle_data:
                column3.AddGate(Gate("I",1,[i]))

            self.circuit.AddColumn(column3)           

            # bottom cnot
            column4 = Column.Column(n+a)
            for i in range(a):
                column4.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][2]]))

            # idle on qubits not touched
            all_data = set(range(n))
            for i in range(len(self.CheckMatrix)):
                all_data.discard(self.CheckMatrix[i][2])
            idle_data = list(all_data)
            for i in idle_data:
                column4.AddGate(Gate("I",1,[i]))

            self.circuit.AddColumn(column4)
            
            
            # left cnot
            column5 = Column.Column(n+a)
            for i in range(a):
                column5.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][3]]))
                
            # idle on qubits not touched
            all_data = set(range(n))
            for i in range(len(self.CheckMatrix)):
                all_data.discard(self.CheckMatrix[i][3])
            idle_data = list(all_data)
            for i in idle_data:
                column5.AddGate(Gate("I",1,[i]))
            self.circuit.AddColumn(column5)
                            

            # measure ancilla
            column6 = Column.Column(n+a)
            for i in range(n,n + a):
                column6.AddGate(Gate("MX",1,[i]))
            #idle on data
            for i in range(n):
                column6.AddGate(Gate("I",1,[i]))
                
            self.circuit.AddColumn(column6)