import sys, os,argparse
from pathlib import Path
from tempfile import NamedTemporaryFile
from math import sqrt

import numpy as np
import parmed as pmd
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

try:
    import openmm
    from openmm import unit
except ImportError:
    from simtk import openmm, unit

import openff.toolkit
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField

from pdbfixer import PDBFixer


def fix_pdb(pdb, base, rns):
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if (rns):
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
    fixer.addMissingHydrogens(7.0)
    openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open(f"{base}_prep.pdb", 'w'))


def parm_complex(prot, lig, ff, wm):
    """[summary]

    Args:
        prot (pdb): [description]
        lig (sdf): [description]
        ff (protein ff): [description]
        wm (water model): [description]
    """
    receptor_path = prot
    ligand_path = lig

    ligand = Molecule.from_file(ligand_path)
    ligand_positions = ligand.conformers[0]
    ligand_topology = ligand.to_topology()

    omm_forcefield = openmm.app.ForceField(ff, wm)
    smirnoff = SMIRNOFFTemplateGenerator(forcefield="openff-2.0.0.offxml", molecules=[ligand])
    omm_forcefield.registerTemplateGenerator(smirnoff.generator)

    pdb = openmm.app.PDBFile(receptor_path)
    modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand_topology.to_openmm(), ligand_positions)

    modeller = solvate(modeller, omm_forcefield)

    system = omm_forcefield.createSystem(modeller.topology, nonbondedMethod=openmm.app.PME, rigidWater=False)
    topology = modeller.getTopology()

    complex_struct = pmd.openmm.load_topology(topology, system, xyz=modeller.getPositions())
    complex_struct.save('complex.prmtop', overwrite=True)
    complex_struct.save('complex.inpcrd', overwrite=True)


def solvate(modeller, ff):
    boxSize, padding, boxVectors = None, None, None
    geompadding = float(6) * unit.nanometer
    maxSize = max(max((pos[i] for pos in modeller.getPositions())) - min((pos[i] for pos in modeller.getPositions())) for i in range(3))
    vectors = openmm.Vec3(1, 0, 0), openmm.Vec3(1/3, 2*sqrt(2)/3, 0), openmm.Vec3(-1/3, sqrt(2)/3, sqrt(6)/3)
    boxVectors = [(maxSize+geompadding)*v for v in vectors]
    modeller.addSolvent(ff, model='tip3p', boxVectors=boxVectors)
    return modeller


def assert_watermodel(wm):
    assert wm in ['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'spce'], 'Unknown water model %s\n' % wm
    wm = 'amber14/%s.xml' % wm
    return wm


def assert_forcefield(ff):
    assert ff in ['ff14SB', 'ff15ipq'], 'Unknown forcefield %s\n' % ff
    ff = 'amber14/protein%s.xml' % ff
    return ff


def assert_extraforcefield(effs):
    effs_w_paths = []
    for ff in effs:
        assert ff in ['DNA.OL15', 'DNA.bsc1', 'RNA.OL3', 'lipid17', 'GLYCAM_06j-1'], 'Unknown extra forcefield %s\n' % ff
        effs_w_paths.append('amer14/%s.xml' % ff)
    return effs_w_paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates pre-production files for AMBER \
    MD. Can handle a receptor or ligand by themselves (checks for ligand library \
    files in the current directory and generates them if they don't exist) or \
    sets up the complex if given both a receptor and ligand.")

    parser.add_argument('-s', '--structures', nargs='+', required=True, help='Structures \
    for which you want to run a simulation. N.B. if more than one is provided \
    they will be simulated together.')

    parser.add_argument('-n', '--out_name', 
    help='Optionally provide a filename prefix for the output.  For multi-structure inputs, default is "complex."')

    parser.add_argument('-p', '--libs', nargs='+', required=False, help="Optionally specify \
    a prefix for nonstandard residue library files; this can include their path if they aren't \
    in the current directory. If the relevant files exist we assume you want \
    to use them, otherwise we assume this is where you want them to go and \
    derive the residue name accordingly.")

    parser.add_argument('-w', '--water_dist', default=12, help='Water box \
    distance; defaults to 12.')

    parser.add_argument('-wm', '--water_model', default='tip3p',
            help='Water model; tip3p, tip3pfb, tip4pew, tip4pfb, and spce are available. \
            Defaults to tip3p.')

    parser.add_argument('-ff', '--force_field', default='ff15ipq',
            help='Force field; ff15ipq and ff14SB are available. \
                Defaults to ff15ipq.')

    parser.add_argument('-eff', '--extra_force_field', nargs='+', default='',
            help='Extra force fields (e.g. DNA, lipids); DNA.OL15, DNA.bsc1, RNA.OL3, lipid17, and GLYCAM_06j-1 are available. \
                defaults to null.')

    parser.add_argument('-t', '--temperature', default=300, help='Simulation \
    temperature; defaults to 300K.')

    parser.add_argument('-l', '--prod_length', default=100, help='Length of \
    production-run MD. This is used to generate the input files for that MD, \
    but note that by default this MD is not run. Units are nanoseconds.')

    parser.add_argument('-k', '--keep_velocities', default=False,
            action='store_true', help='Keep velocities from preproduction run \
                    when starting the production MD. Default is false.')

    parser.add_argument('-r', '--run_prod_md', default=False,
    action='store_true', help='Run production MD locally when all preprocessing \
    is finished. Default is false, because you might want to run it on a \
    cluster.')

    parser.add_argument('-noh', '--no_touch_hyd', dest='noh', default=False,
            action='store_true', help="Don't remove any hydrogens.")

    parser.add_argument('-nc', '--net_charge', help='Optionally specify a net \
            charge for small molecule parametrization with antechamber.')

    parser.add_argument('-parm', '--parm_only', action='store_true', default =
    False, help="Only generate the necessary ligand parameters, don't do the \
    preproduction MDs")

    parser.add_argument('-ui', '--uninteractive', action='store_true',
            default=False, help="Turn off interactive mode, which exists to let \
            you check the output of pdb4amber for potentially serious problems \
            with input structures.")

    parser.add_argument('-O', '--overwrite', action='store_false',
            default=True, help='Overwrite files in the current directory that \
            are generated by this script. Default is True.')

    parser.add_argument('-df', '--coord_dump_freq', default=5000,
    help='Frequency for dumping coordinates to traj_file. Defaults to 5000. \
    The old script referred to this as the "timestep."')
    
    args = parser.parse_args()

    # Assert water model and forcefields are in openmm
    args.water_model = assert_watermodel(args.water_model)
    args.force_field = assert_forcefield(args.force_field)
    args.extra_force_field = assert_extraforcefield(args.extra_force_field)
    
    for structure in args.structures:
        assert os.path.isfile(structure),'%s does not exist\n' % structure
        struct_file_tup = os.path.splitext(structure)
        base = struct_file_tup[0]
        ext = struct_file_tup[1] 
        # fix_pdb(structure, base, True)
        parm_complex(f"{base}_prep.pdb", 'gfp_ligand.sdf', args.force_field, args.water_model)


            

        

    
