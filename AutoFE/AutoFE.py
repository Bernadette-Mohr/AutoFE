import argparse
from pathlib import Path
import graph2topology
import graph2coordinates
import minimize_graphs
import equilibrate_charge
import copy_files
import system_setup
import windows


def main():
    """
     mols_archive: pathlib Path to pickle output of networkx graphs to process
     environment_path: pathlib Path of the Gromacs input files directory for the environments
     mdp_path: pathlib Path to directory containing the mdp options files
     forcefield_file: pathlib Path to the forcefield file used for generating the small molecules
    """
    parser = argparse.ArgumentParser(description='Automatic free energy calculations of arbitrary small molecules.')
    parser.add_argument('-mol', '--molecules', type=Path, required=True,
                        help='Need: location of molecule graphs archive')
    parser.add_argument('-env', '--environment', type=Path, required=True,
                        help='Need: location of the Gromacs input files for the simulation environments.')
    parser.add_argument('-mdp', '--parameters', type=Path, required=True,
                        help='Location of mdp options files. Prepend \"VACUUM-\" for small molecule energy '
                             'minimization and \"WATER-\" for simulation in water box.')
    parser.add_argument('-ff', '--forcefield', type=Path, required=True,
                        help='Location of main forcefield file to be used for molecule generation.')
    parser.add_argument('-c', '--cluster', type=Path, required=True,
                        help='Need: location of setup and run control files for running simulations on a cluster.')

    args = parser.parse_args()

    mols_archive = args.molecules
    environment_path = args.environment
    mdp_path = args.parameters
    forcefield_file = args.forcefield
    run_path = args.cluster

    """
     Create Gromacs topologies for the molecule graphs found in the pickle archive. Returns parent directory of the 
     generated molecule directories.
    """
    molecules_path, molecules_info = graph2topology.create_topology(mols_archive, forcefield_file)

    """
     Place all particles of the molecule randomly inside the volume of a sphere around the box center, make sure they 
     have at least a distance of the shortest bond length to each other.
    """
    graph2coordinates.create_coordinates(molecules_path, molecules_info)

    """
     Energy minimization of the small molecules in vacuum to restore the structure from the bonds and angles in the 
     topology. Minimization routine adapted from "backward.py" 
     (http://www.cgmartini.nl/index.php/downloads/tools/240-backward).
    """
    minimize_graphs.run_minimization(molecules_path, mdp_path)

    """
     Adjust the type/number of charges present in the environment to equilibrate any charge present in the small 
     molecule. Returns the paths for the adapted environment.gro files and information about presence of ions in the 
     adapted environment.
    """
    env_paths, mol_charge, contains_ions = equilibrate_charge.equilibrate(environment_path, molecules_path, molecules_info)

    """
     Copy MD parameter files and molecule topology and coordinate files to the respective simulation directories. 
     Insert the name of the environment and the correct coupling groups in the mdp files.
    """
    copy_files.copy_files(env_paths, molecules_path, args.parameters, run_path, mol_charge=mol_charge, contains_ions=contains_ions)

    """
     Create the initial system coordinate, index and topology files for simulating each molecule in every environment.
    """
    system_setup.setup_system(env_paths)

    """
     Set up subdirectories for simulating all lambda steps for each molecule and environment.
    """
    windows.create(env_paths, forcefield_file, mol_charge)

    """
     Submits all lambda step subdirectories in all environment simulation directories pairwise to a cluster. The path 
     to the cluster environment and cluster run control files have to be provided via the -c/--cluster flag.
    """


if __name__ == "__main__":
    main()
