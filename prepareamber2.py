import argparse
from pdbfixer import PDBFixer
import os, glob


'''
Strip path and extension from filename in order to generate other filenames
'''
def get_base(fname):
    return os.path.splitext(os.path.basename(fname))[0]



def fix_pdb(structures):
    for structure in structures:
        if 'pdb' in os.path.splitext(structure)[-1]:
            # the "structure" string in the args.structure list will be updated so
            # that it corresponds to the PDB we should use for subsequent steps
            assert os.path.isfile(structure),'%s does not exist\n' % structure

            # use pdbfixer to find anything wrong with initial pdb(s)
            fixer = PDBFixer(filename=structure)
            # find missing residues, atoms
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            # add missing residues and atoms
            print('Adding missing residues and atoms...')
            fixer.addMissingAtoms()

            # TODO: deal with nonstandard residues
            fixer.findNonstandardResidues()

            # if modifications files not provided
            # if not args.libs:
            #     args.libs = []
            # find_nonstandard_libraries(structure, args.libs, fixer.nonstandardResidues)
            # for res in fixer.nonstandardResidues:
            #     print('Unknown residue {name} at {id} in {chain}. '\
            #         'Suggested replacement: {replacement}'.format(
            #         name=res[0].name, id=res[0].id, chain=res[0].chain, replacement=res[1]))

            


def main():
    parser = argparse.ArgumentParser(description="Generates pre-production files for AMBER \
    MD. Can handle a receptor or ligand by themselves (checks for ligand library \
    files in the current directory and generates them if they don't exist) or \
    sets up the complex if given both a receptor and ligand.")

    parser.add_argument('-s', '--structures', nargs='+', required=True, help='Structures \
    for which you want to run a simulation. N.B. if more than one is provided \
    they will be simulated together.')

    parser.add_argument('-p', '--libs', nargs='+', required=False, help="Optionally specify \
    a prefix for nonstandard residue library files; this can include their path if they aren't \
    in the current directory. If the relevant files exist we assume you want \
    to use them, otherwise we assume this is where you want them to go and \
    derive the residue name accordingly.")

    parset.add_argument('-ns', '--replace_nonstandard', default=False, action='store_true',
    help="Replace nonstandard residues with suggested replacement if nonstandard residue libraries \
    are not found and not provided. Deafult is false, becuase you may not have known there are \
    missing libraries.")

    args = parser.parse_args()

    fix_pdb(args.structures)

if __name__ == '__main__':
    main()