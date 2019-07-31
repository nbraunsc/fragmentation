import sys
from sys import argv
import Avogadro as av

mol = av.MoleculeFile.readMolecule('carbonylavo.cml')

for bond in mol.bonds:
    first = bond.beginAtom
    last = bond.endAtom
    print(first.index, last.indec, bond.length)
