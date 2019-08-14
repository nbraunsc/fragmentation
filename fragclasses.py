import os
import copy
import numpy as np
import sys
from sys import argv
import xml.etree.ElementTree as ET

class Molecule():
    def __init__(self):
        #number of atoms
        self.natoms = {} 
        #atom coords
        self.atomsxyz = []
        #atom element type
        self.atomselement = []
        #atom indexs
        self.atomsindex = []
        #atom coords with index and element, type:list
        self.atomtable = []
        #table with labels of integers for atoms, starting at 0, type: np array
        self.bond_lables = []
        #table with just bond orders, type: np array
        self.bond_order = []
        #table with labels and bond order, type:numpy.array
        self.bond_table = []
        #connectivity matrix
        self.A = []
        #single bond list
        self.sbond = []
        #double bond list
        self.dbond = []
        #list of lists for primiatives
        self.prims = []
        #first set of frags
        self.frags = []
        #list of unique fragments:
        self.uniquefrags = []
        self.frag_inter = []

    def parse_cml(self, filename):
        self.filename = filename
        tree = ET.parse(filename)
        self.tree = tree
        root = tree.getroot()
        molecule = root
        self.molecule = molecule
        atomArray = root[0]
        self.atomArray = atomArray
        bondArray = root[1]
        self.bondArray = bondArray        
        self.natoms = len(self.atomArray)
        self.A = np.zeros( (self.natoms,self.natoms)) 
        #start of pulling molecule info from cml file
        self.atomtable = [[0 for x in range(4)] for y in range(self.natoms)]
        for atomi in range(0, self.natoms):
            self.atomsxyz = list([float(atomArray[atomi].attrib['x3']) , float(atomArray[atomi].attrib['y3']) , float(atomArray[atomi].attrib['z3'])])
            self.atomselement = atomArray[atomi].attrib['elementType']
            #self.atomsindex = int(atomi)
            self.atomtable[atomi][0] = self.atomselement
            #self.atomtable[atomi][1] = self.atomsindex
            self.atomtable[atomi][1] = float(atomArray[atomi].attrib['x3'])
            self.atomtable[atomi][2] = float(atomArray[atomi].attrib['y3'])
            self.atomtable[atomi][3] = float(atomArray[atomi].attrib['z3'])

        #start extracting bond order
        for bondi in bondArray:
            a12 = bondi.attrib['atomRefs2'].split()
            #assert is making all first part of indexes = a
            assert(a12[0][0] == "a")
            assert(a12[1][0] == "a")
            #takes a1 and makes it a0 so index of element is same in bondorder
            a1 = np.array([int(a12[0][1:])-1])
            a2 = np.array([int(a12[1][1:])-1])
            self.bond_labels = np.append(a1,a2)
            self.bond_order = np.array([int(bondi.attrib["order"])])
            self.bond_table = np.append(self.bond_labels, self.bond_order)
            #these need to be within for loop, don't know if need self.
            x = self.bond_table[0]
            y = self.bond_table[1]
            z = self.bond_table[2]
            self.A[x][y] = z
            self.A[y][x] = z
            
        #returning indices in connectivity chart where bond order is 1 or 2 in form of (row, colum)
        single = np.where(self.A == 1)
        double = np.where(self.A > 1)
        self.sbond = list(zip(single[0], single[1]))
        self.dbond = list(zip(double[0], double[1]))

    #makes list of list with primatives, have some repeating primiatives
    def get_prims(self):
        for i in range(0, len(self.atomtable)):
            if self.atomtable[i][0] != "H":
                self.prims.append([i])
                for j in range(0, len(self.atomtable)):
                    if self.A[i][j] == 1 and self.atomtable[j][0] == "H" or self.A[i][j] > 1:
                        self.prims[-1].append(j)
                    for k in range(0, len(self.atomtable)):
                        if self.A[i][j] > 1 and self.A[j][k] == 1 and self.atomtable[k][0] != "C":
                            self.prims[-1].append(j)
                            self.prims[-1].append(k)
        
        #deletes duplicates within a prim
        self.prims = list(set(x) for x in self.prims)

        #gets rid of repeating prims, sorted makes sure the order doesn't matter
        for i in range(0, len(self.prims)):
            self.prims[i] = tuple(sorted(self.prims[i]))
        self.prims = set(self.prims)
        
    #eta is how high you want to go, eta = 0: prims, eta = 1:one bond, eta = 2:two bonds etc., iteration = 0 to start
    def get_frags(self, eta, iteration):
        self.newfrag =[]
        if iteration == 0:
            for prim in self.prims:
                self.frags.append(list(prim))
                
        if iteration != 0:
            for frag in self.frags:
                self.newfrag.append([])
                for atom in frag:
                    for prim in self.prims:
                        for atom2 in prim:
                            if self.A[atom][atom2] != 0:
                                if atom2 in frag:
                                    continue
                                for atom3 in prim:
                                    self.newfrag[-1].append(atom3)
        for i in range(0, len(self.newfrag)):
            for atom in self.newfrag[i]:
                self.frags[i].append(atom)
        if iteration < eta:
            iteration = iteration + 1
            self.get_frags(eta, iteration)

    def remove_frags(self):
        self.uniquefrags = []
        self.frags = list(set(x) for x in self.frags)
        print(self.frags)
        for i in range(0, len(self.frags)):
            
            for j in range(i+1, len(self.frags)):
                print(str(i) + " " +str(j)) 
                if self.frags[i].issubset(self.frags[j]):
                    continue
                self.uniquefrags.append(self.frags[i])
                    #if all atoms are described with eta bonds, remove all the remaining fragments
                if self.natoms == len(self.uniquefrags[0]):
                    for k in range(1, len(self.uniquefrags)):
                        self.uniquefrags.remove(self.uniquefrags[k])
        #if eta is larger than bonds allow thus all i's are a subset of j add 1st entry in self.frags
        if len(self.uniquefrags) == 0:
            self.uniquefrags.append(self.frags[0])
    
    def frag_conn(self):
        self.frag_inter = []
        for i in self.uniquefrags:
            for j in self.uniquefrags:
                if i != j:
                    if i.isdisjoint(j):
                        continue
                    self.frag_inter = i.intersection(j)

if __name__ == "__main__":
    carbonyl = Molecule()
    carbonyl.parse_cml("largermol.cml")
    carbonyl.get_prims()
    carbonyl.get_frags(1, 0)
    carbonyl.remove_frags()
    carbonyl.frag_conn()
    print(carbonyl.uniquefrags)
