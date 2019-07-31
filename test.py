import os
import numpy
import sys
from sys import argv
import xml.etree.ElementTree as ET

class Molecule:
    def __init__(self):
        #number of atoms
        self.natoms = 0
        #atom objects/element names
        self.atoms = []
        #atom coords.xyz
        self.atomcoord = []
        self.bond_table = []
        self.bond_order = {}
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
            self.atoms.append([])
            self.atoms[-1].element = atomArray[atomi].attrib['elementType']
            self.atoms[-1].index = atomi
            self.bond_table.append([])
            self.atoms[-1].xyz = numpy.asarray([ float(atomArray[atomi].attrib['x3']) , float(atomArray[atomi].attrib['y3']) , float(atomArray[atomi].attrib['z3'])])





carbonyl = Molecule()
carbonyl.parse_cml("carbonylavo.cml")
