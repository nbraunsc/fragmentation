import os
import numpy as np
from sys import argv
import xml.etree.ElementTree as ET

tree = ET.parse("carbonylavo.cml")
root = tree.getroot()
molecule = root
atomArray = root[0]
bondArray = root[1]
natoms = len(atomArray)

atomtable = [[0 for x in range(5)] for y in range(natoms)] 
for atomi in range(0, natoms):
    atomsxyz = list([float(atomArray[atomi].attrib['x3']) , float(atomArray[atomi].attrib['y3']) , float(atomArray[atomi].attrib['z3'])])
    atomselement = atomArray[atomi].attrib['elementType']
    atomsindex = int(atomi)
    atomtable[atomi][0] = atomselement 
    atomtable[atomi][1] = atomsindex 
    atomtable[atomi][2] = float(atomArray[atomi].attrib['x3'])
    atomtable[atomi][3] = float(atomArray[atomi].attrib['y3'])
    atomtable[atomi][4] = float(atomArray[atomi].attrib['z3'])

print(atomtable[0])

A = np.zeros( (natoms,natoms))
for bondi in bondArray:
    a12 = bondi.attrib['atomRefs2'].split()
    assert(a12[0][0] == "a")
    assert(a12[1][0] == "a")
    a1 = np.array([int(a12[0][1:])-1])
    a2 = np.array([int(a12[1][1:])-1])
    bond_labels = np.append(a1,a2)
    bond_dict = np.array([int(bondi.attrib["order"])])
    bond_table = np.append(bond_labels, bond_dict)
    #making connectivity matrix
    x = bond_table[0]
    y = bond_table[1]
    z = bond_table[2]
    A[x][y] = z

#returning indices that have 1 or 2 in their spot (row, column)
single = np.where(A == 1)
double = np.where(A > 1)
listofdouble = list(zip(double[0], double[1]))
listofsingle = list(zip(single[0], single[1]))

#NEED to figure out how to remove duplicates!!!!
for sbond in listofsingle:
    print(type(sbond))      
for dbond in listofdouble:
    print(dbond)             



