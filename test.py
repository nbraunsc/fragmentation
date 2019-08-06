import os
from sys import argv
from fragclasses import *

carbonyl = Molecule()
carbonyl.parse_cml("carbonylavo.cml")
carbonyl.get_prims()
print(carbonyl.prims)
carbonyl.get_frags(3, 0)
print(carbonyl.frags)
