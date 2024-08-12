# Author: Diego Ruiz
# Affiliation: Alice&Bob - INRIA
# Date: 2023

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
import matplotlib
import numpy as np
import qiskit as qt
import random
import Column
import matplotlib.pyplot as plt
import timeit
from Gate import Gate
from numpy.random import *
import copy

# class to represent the quantum circuit that contains columns
class Circuit:
    
    def __init__(self,nqubit):
        # list of columns composing the circuit
        self.columns = []
        # number of qubits
        self.n = nqubit
        # number of columns
        self.nbcol = 0
        # current errors on the qubits 
        self.psi = [0]*nqubit
        # store the psi at all steps of the circuit
        self.AllPsi = []
        # list of measurements 
        self.Mmt = []
              
    # add column to the circuit
    def AddColumn(self,column):
        # check column matches circuit
        if column.n != self.n :
            print("Column has not the same number of qubits as the circuit")
            
        self.columns.append(column)
        self.nbcol+=1

    # Trigger errors in the circuit and propagate them
    def Error(self,bprint):
        for i in range(self.nbcol):
            column = self.columns[i]
            # propagate errors from previous step
            column.Propagate(self.psi,self.Mmt)
            # trigger errors
            column.Error(self.psi) 
            # store current errors
            if bprint:
                self.AllPsi.append(copy.deepcopy(self.psi))

    # Reset the circuit for a new run            
    def Reset(self):
        self.Mmt = []
        self.AllPsi = []
        self.psi = [0]*self.n

    # print the circuit with qiskit
    def print(self,k,fileloc):
                
        CounterMmt = 0
        
        # Build a quantum circuit
        Qiskitcircuit = QuantumCircuit(self.n, 100)
        # Custom gate
        gateprep = qt.circuit.Gate("P", 1, [])
        gatemeasure = qt.circuit.Gate("M",1,[])
        gateIdle = qt.circuit.Gate("I", 1, [])
        
        for j in range(self.nbcol):
            
            for i in range(self.columns[j].nbgates):
                
                gate = self.columns[j].gates[i]
                
                # Qiskit gate
                if gate.name == "Z" :
                    Qiskitcircuit.z(gate.qubits[0])
                if gate.name == "X" :
                    Qiskitcircuit.x(gate.qubits[0])
                if gate.name == "CX" : 
                    Qiskitcircuit.cx(gate.qubits[0],gate.qubits[1])
                if gate.name == "CCX" : 
                    Qiskitcircuit.ccx(gate.qubits[0],gate.qubits[1],gate.qubits[2])
                if gate.name == "CZ" : 
                    Qiskitcircuit.cz(gate.qubits[0],gate.qubits[1])
                
                # Custom gate
                if gate.name == "P" : 
                    Qiskitcircuit.append(gateprep,gate.qubits)
                if gate.name == "I" : 
                    Qiskitcircuit.append(gateIdle,gate.qubits)
                if gate.name == "MX" :
                    Qiskitcircuit.append(gatemeasure,gate.qubits)
                    if len(self.Mmt)>0 and self.Mmt[CounterMmt]==1:
                        Qiskitcircuit.z(gate.qubits[0])
                    CounterMmt+=1
                        
            # barrier between columns           
            Qiskitcircuit.barrier(range(self.n))
            
            # print psi
            if len(self.AllPsi)>j:
                for i in range(len(self.AllPsi[j])):
                    if self.AllPsi[j][i]==1:
                        Qiskitcircuit.z(i)
                        
                Qiskitcircuit.barrier(range(self.n))
              
        # draw the circuit
        Qiskitcircuit.draw(output='mpl',
                           style={'displaycolor': {'M': ('#000000', '#FFFFFF'),
                                                   'z':('#F00000', '#000000'),
                                                   'I':('#D3D3D3', '#000000'),
                                                   'ccx':('#000000', '#000000')}},
                           scale=1,
                           filename=fileloc+str(k)+'.png',
                           justify='right',fold=300,
                           cregbundle=True)
        
        plt.close()

        
