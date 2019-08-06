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
        self.bond_dict = []
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
        self.atomtable = [[0 for x in range(5)] for y in range(self.natoms)]
        for atomi in range(0, self.natoms):
            self.atomsxyz = list([float(atomArray[atomi].attrib['x3']) , float(atomArray[atomi].attrib['y3']) , float(atomArray[atomi].attrib['z3'])])
            self.atomselement = atomArray[atomi].attrib['elementType']
            self.atomsindex = int(atomi)
            self.atomtable[atomi][0] = self.atomselement
            self.atomtable[atomi][1] = self.atomsindex
            self.atomtable[atomi][2] = float(atomArray[atomi].attrib['x3'])
            self.atomtable[atomi][3] = float(atomArray[atomi].attrib['y3'])
            self.atomtable[atomi][4] = float(atomArray[atomi].attrib['z3'])

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
            self.bond_dict = np.array([int(bondi.attrib["order"])])
            self.bond_table = np.append(self.bond_labels, self.bond_dict)
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
                    if self.A[i][j] != 0 and self.atomtable[j][0] == "H" or self.A[i][j] > 1: 
                        self.prims[-1].append(j)
        #gets rid of repeating prims, sets make sure order doesn't matter
        for i in range(0, len(self.prims)):
            self.prims[i] = tuple(sorted(self.prims[i]))
        self.prims = set(self.prims)

    #eta is how high you want to go order wise, iteration = 0 to start
    def get_frags(self, eta, iteration):
        for prim in self.prims:
            for frag in self.frags:
                if prim in frag:
                    continue
                self.frags.append(list(prim))
            for atom in prim:
                for prim2 in self.prims:
                    for atom2 in prim2:
                        if self.A[atom][atom2] != 0:
                            if atom2 in self.frags:
                                continue
                            for atom3 in prim2:
                                self.frags.append(prim2)
                                #self.frags[-1].append(prim2)
        #for i in range(0, len(self.frags)):
            #self.frags[i] = set(self.frags[i])
        iteration = iteration + 1
        if iteration < eta:
            print('d')
            self.get_frags(eta, iteration)

if __name__ == "__main__":
    carbonyl = Molecule()
    carbonyl.parse_cml("carbonylavo.cml")
    carbonyl.get_prims()
    carbonyl.get_frags(1, 0)
    print(carbonyl.frags)


