# from openmm.app import *
# from openmm import *
# from openmm.unit import *
from pdbfixer import *


def main():
    fixer = PDBFixer(filename='complex_noH.pdb')
    fixer.findNonstandardResidues()
    print(fixer.nonstandardResidues)


if __name__ == '__main__':
    main()