import os
import numpy as np
import sys
from sys import argv
import xml.etree.ElementTree as ET

class Molecule:
    def __init__(self):
        #number of atoms
        self.natoms = 0
        #atom coords
        self.atomsxyz = []
        #atom element type
        self.atomselement = []
        #atom indexs
        self.atomsindex = []
        #numpy array of element type and index
        self.atomselement_index = []
        #atom coords with index and element, type: numpyarray
        self.atomtable = []
        #table with labels of integers for atoms, starting at 0, type:numpy.array
        self.bond_lables = []
        #table with just bond orders, type: numpy.array
        self.bond_dict
        #table with labels and bond order, type:numpy.array
        self.bond_table = []
        Atom = open('atomcoords.txt','a+')

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
        self.natoms = len(atomArray)
        
        #start of pulling molecule info from cml file
        for atomi in range(0, self.natoms):
            self.atomsxyz = np.asarray([ float(atomArray[atomi].attrib['x3']) , float(atomArray[atomi].attrib['y3']) , float(atomArray[atomi].attrib['z3'])])
            self.atomselement = np.array([atomArray[atomi].attrib['elementType']])
            self.atomsindex = np.array([int(atomi)])
            #numpy.append is adding element to index in a column
            self.atomselement_index = np.append(self.atomselement, self.atomsindexes)
            #numpy.append is adding element-index to coordinates of molecule in column
            self.atomtable = np.append(self.atomselement_index, self.atomsxyz)
#A needs to be defined outside of for loop and class??
A = numpy.zeros( (natoms,natoms))
        
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
            self.bond_table = list(np.append(self.bond_labels, self.bond_dict))
            #these need to be within for loop, don't know if need self.
            x = bond_table[0]
            y = bond_table[1]
            z = bond_table[2]
            A[x][y] = z
            A[y][x] = z
#returning indices in connectivity chart where bond order is 1 or 2 in form of (row, colum)
single = np.where(A == 1)
double = np.where(A > 1)
listofdouble = list(zip(double[0], double[1]))
listofsingle = list(zip(single[0], single[1]))
for sbond in listofsingle:
    singlebond = np.array([sbond])
    print(sbond)
for dbond in listofdouble:
    doublebond = np.array([dbond])
    print(dbond)


carbonyl = Molecule()
carbonyl.parse_cml("carbonylavo.cml")




