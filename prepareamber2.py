import argparse
from pdbfixer import PDBFixer
import pdb_util as util
import os


def fix_pdb(structures):
    for strucutre in structures:
        if 'pdb' in os.path.splitext(structure)[-1]:
            #the "structure" string in the args.structure list will be updated so
            #that it corresponds to the PDB we should use for subsequent steps
            assert os.path.isfile(pdb),'%s does not exist\n' % structure

            #use pdbfixer to find anything wrong with initial pdb(s)
            fixer = PDBFixer(filename=structure)
            fixer.findMissingResidues()
            print(fixer.missingResidues)
            fixer.addMissingAtoms()
        

def main():
    parser = argparse.ArgumentParser(description="Generates pre-production files for AMBER \
    MD. Can handle a receptor or ligand by themselves (checks for ligand library \
    files in the current directory and generates them if they don't exist) or \
    sets up the complex if given both a receptor and ligand.")

    parser.add_argument('-s', '--structures', nargs='+', required=True, help='Structures \
    for which you want to run a simulation. N.B. if more than one is provided \
    they will be simulated together.')

    args = parser.parse_args()

    fix_pdb(args.structures)

if __name__ == '__main__':
    main()