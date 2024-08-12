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
              
    # Parity check matrix of the code
    def CreateCheckMatrix(self):
       
        self.CheckMatrix = np.zeros(((self.H-1)*self.L,4),dtype = int )
        
        L = self.L
       
        for i in range(self.H-1):
            for j in range(self.L):
                self.CheckMatrix[i*self.L+j] = np.array([(i+1)*L+j,
                                                         i*L + (j+1)%self.L,
                                                         i*L + j,
                                                         i*L + (j-1+self.L)%self.L])
                
    # create the circuit of the cellular automaton code from the parity check matrix
    def CreateCircuit(self,d):
        
        n = self.H*self.L # nb of data
        a = (self.H-1)*self.L # nb of ancilla
        
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
            self.circuit.AddColumn(column2)

            
            # right cnot
            column3 = Column.Column(n+a)
            for i in range(a):
                column3.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][1]])) 
            self.circuit.AddColumn(column3)

            
            # bottom cnot
            column4 = Column.Column(n+a)
            for i in range(a):
                column4.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][2]]))
            self.circuit.AddColumn(column4)


            # left cnot
            column5 = Column.Column(n+a)
            for i in range(a):
                column5.AddGate(Gate("CX",2,[n+i,self.CheckMatrix[i][3]]))

                
            self.circuit.AddColumn(column5)
            
            # measure ancilla
            column6 = Column.Column(n+a)
            for i in range(n,n + a):
                column6.AddGate(Gate("MX",1,[i]))
                
            self.circuit.AddColumn(column6)
                
        
        

        
