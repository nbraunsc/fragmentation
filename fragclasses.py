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
        self.primchart = []
        self.molchart = []
        #first set of frags
        self.fragments = []
        #unqiue fragments
        self.uniquefrags = []
        self.fullmolecule = []

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
                        if self.A[i][j] > 1 and self.A[j][k] == 1 and self.atomtable[k][0] == "H":
                            self.prims[-1].append(j)
                            self.prims[-1].append(k)
        
        #deletes duplicates within a prim
        self.prims = list(set(x) for x in self.prims)
        for i in range(0, len(self.prims)):
            self.prims[i] = tuple(sorted(self.prims[i]))
        self.prims = set(self.prims)
    
    #if spot in A is non zero add 1 to primchart in row and column of prim1 and prim2 
    def prim_conn(self):
        self.prims = list(self.prims)
        self.primsleng = len(self.prims)
        self.primchart = np.zeros( (self.primsleng,self.primsleng))
        for prim1 in range(0, len(self.prims)):
            for atomi in self.prims[prim1]:
                for prim2 in range(0, len(self.prims)):
                    for atomj in self.prims[prim2]:
                        if prim1 == prim2:
                            continue
                        if self.A[atomi][atomj] != 0:
                            self.primchart[prim1][prim2] = 1
                            self.primchart[prim2][prim1] = 1
    
    def get_molmatrix(self, eta, i):    #i must 2 to start
        if i == 2:
            self.molchart = self.primchart.dot(self.primchart)
            np.fill_diagonal(self.molchart, 0)
            for x in range(0, len(self.prims)):
                for y in range(0, len(self.prims)):
                    if self.molchart[x][y] != 0:
                        self.molchart[x][y] = i
            
            for x in range(0, len(self.prims)):      #checkpoint, unsure if needed
                for y in range(0, len(self.prims)):
                    if self.primchart[x][y] != 0:
                        self.molchart[x][y] = 0
            self.molchart = np.add(self.molchart, self.primchart)
            
        if i > 2:
            self.molchartnew = self.molchart.dot(self.primchart)
            np.fill_diagonal(self.molchartnew, 0)
            for x in range(0, len(self.prims)):
                for y in range(0, len(self.prims)):
                    if self.molchartnew[x][y] != 0:
                        self.molchartnew[x][y] = i
            
            for x in range(0, len(self.prims)):     #checkpoint against previous matrix
                for y in range(0, len(self.prims)):
                    if self.molchart[x][y] != 0:
                        self.molchartnew[x][y] = 0
            self.molchart = np.add(self.molchart, self.molchartnew)

        if i < eta:     #recursive part of function
            i = i+1
            self.get_molmatrix(eta, i)
    
    def get_frag(self, deg):    #deg is the degree of fragments wanted
        self.frag = []
        for x in range(0, len(self.molchart)):
            for y in range(0, len(self.molchart)):
                if self.molchart[x][y] <= deg and self.molchart[x][y] != 0:
                    self.frag.append([x])
                    self.frag[-1].append(y)
        
        for z in range(0, len(self.frag)):
            for w in range(0, len(self.frag)):
                if z == w:
                    continue
                if self.frag[z][0] == self.frag[w][0]:
                    self.frag[z].extend(self.frag[w][:])    #combines all prims with frag connectivity <= eta, list of lists
        
        for i in range(0, len(self.frag)): 
            self.frag[i] = set(self.frag[i])    #makes into a list of sets
        
        for k in range(0, len(self.frag)):
            if self.frag[k] in self.fragments:
                continue
            self.fragments.append(self.frag[k])
                   
    def remove_frags(self):
        for i in range(0, len(self.fragments)):
            add = True
            for j in range(i+1, len(self.fragments)):
                if self.fragments[i].issubset(self.fragments[j]) or self.fragments[j].issubset(self.fragments[i]) or self.fragments[i] == self.fragments[j]:
                    add = False
            if add == True:
                self.uniquefrags.append(self.fragments[i])
        
        #for frag in self.fragments:       #in here in case full molecule is described
            #if len(frag) == len(self.prims):
                #self.uniquefrags.append(frag)

    def frag_conn(self):        #finds union of frags, problem is idk how to properly map fragments to correct interstions for calculations
        self.fragconn = []
        for i in range(0, len(self.uniquefrags)):
            add = True
            for j in range(i+1, len(self.uniquefrags)):
                if self.uniquefrags[i].isdisjoint(self.uniquefrags[j]):
                    add = False
                if add == True:
                    self.fragconn.append(self.uniquefrags[i].intersection(self.uniquefrags[j]))

        print(self.fragconn)

    def make_fragxyz(self):
        atomlist = []
        for frag in range(0, len(self.uniquefrags)):
            self.uniquefrags[frag] = list(self.uniquefrags[frag])

        for i in range(0, len(self.uniquefrags)):   #makes new list for each fragment
            self.fragxyzi = []

        for frag in self.uniquefrags:   #adding coordinates to list of lists
            for prim in frag:
                for i in self.prims:
                    for atom in i:
                        atomlist[frag][-1].append([self.prims[prim]])
        print(atomlist)
                
                #self.fragxyzi.append(self.atomtable[j])
    #def make_connxyz(self):


if __name__ == "__main__":
    carbonyl = Molecule()
    carbonyl.parse_cml("aspirin.cml")
    carbonyl.get_prims()
    print(carbonyl.atomtable)
    print(carbonyl.prims)
    carbonyl.prim_conn()
    print(carbonyl.primchart)
    carbonyl.get_molmatrix(10, 3)
    print(carbonyl.molchart)
    carbonyl.get_frag(1)
    carbonyl.remove_frags()
    print(carbonyl.uniquefrags)
    carbonyl.frag_conn()

