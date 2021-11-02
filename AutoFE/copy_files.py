import sys
import shutil
import fileinput
from contextlib import closing
from pathlib import Path
import regex as re
# from typing import Any
from equilibrate_charge import count_ions

"""
 Copy the minimized molecule coordinates and corresponding topologies to their subdirs in the environment directories. 
 Removes the minimization step number from the coordinate file.
"""


def copy_mols(mol_path, new_dir, mol):
    #
    itp_file = Path(mol_path, mol + '.itp')
    top_file = Path(mol_path, mol + '.top')
    gro_file = Path(mol_path, mol + '-4.gro')
    old_gro = Path(new_dir, mol + '-4.gro')
    new_gro = Path(new_dir, mol + '.gro')

    files = [itp_file, top_file, gro_file]

    for f in files:
        try:
            shutil.copy(f, new_dir)
        except IOError as e:
            print('Copying molecule files failed: %s' % e)

    try:
        old_gro.rename(new_gro)
    except OSError:
        print('GRO-file renaming failed')


"""
 Replace the required keywords in the sbatch file: environment name, molecule number, number of lambda steps
"""

def adjust_controls(destination, mol_name, charge):

    if charge == 0:
        steps = '40'
    else:
        steps = '80'

    filepath = destination.parent.resolve()
    workdir = str(filepath.relative_to(Path.home().resolve()))

    with closing(fileinput.FileInput(destination, inplace=True)) as sbatch_file:

        for line in sbatch_file:
            rep = {'WORK_DIR': str(workdir), 'NR': mol_name.split('_')[1], 'N_STEPS': steps}
            rep = dict((re.escape(k), v) for k, v in rep.items())
            pattern = re.compile('|'.join(rep.keys()))
            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))


"""
 Copy files that handle setting of cluster variables and running simulations on cluster.
"""


def copy_controls(cluster_path, new_dir, charge, mol_name, env):
    #
    for f in cluster_path.iterdir():

        if charge != 0:
            if f.is_file() and 'CHARGED' in f.name:
                try:
                    name = f'alchemical-{mol_name}.sbatch' # f.name.split('_')[0] + '-' + mol_name + ".sbatch"
                    destination = new_dir.parents[1] / name
                    destination.write_text(f.read_text())
                    adjust_controls(destination, mol_name, charge)
                except IOError as e:
                    print(f'Copying run control files failed: {e}')
        else:
            if f.is_file() and 'NEUTRAL' in f.name:
                try:
                    name = f'alchemical-{mol_name}.sbatch' # f.name.split('_')[0] + '-' + mol_name + ".sbatch"
                    destination = new_dir.parents[1] / name
                    destination.write_text(f.read_text())
                    adjust_controls(destination, mol_name, charge)
                except IOError as e:
                    print(f'Copying run control files failed: {e}')


"""
 Handles insertion of actual environment names instead of 
 placeholders, sets coupling groups according to the presence/absence 
 of ions. 
"""


def change_mdp(filename, env, ions):
    #
    env_name = env.stem

    """
     If module is run standalone, the environment.gro files are searched for the presence of ions.
    """

    if isinstance(ions, Path):
        coordinates = ions.read_text()
        if 'PNA' in coordinates or 'PCL' in coordinates:
            content = True
        else:
            content = False
    else:
        content = ions

    if env_name != 'OCTANE':

        with closing(fileinput.FileInput(filename, inplace=True)) as mdp_file:

            if env_name == 'WATER':
                for line in mdp_file:
                    sys.stdout.write(line.replace('ENV', 'PW'))
            else:
                for line in mdp_file:

                    if env_name == 'CDL2':
                        if content != 0:
                            rep = {'ENV': '{:5.5}'.format(env_name), 'LOW': '1.5', 'HIGH': '2.5'}
                            rep = dict((re.escape(k), v) for k, v in rep.items())
                            pattern = re.compile('|'.join(rep.keys()))
                            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))
                        else:
                            rep = {'ENV': '{:5.5}'.format(env_name), 'PW_ION': 'PW', 'LOW': '1.5', 'HIGH': '2.5'}
                            rep = dict((re.escape(k), v) for k, v in rep.items())
                            pattern = re.compile('|'.join(rep.keys()))
                            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))

                    else:
                        if content != 0:
                            rep = {'ENV': '{:5.5}'.format(env_name), 'LOW': '1.0', 'HIGH': '2.0'}
                            rep = dict((re.escape(k), v) for k, v in rep.items())
                            pattern = re.compile('|'.join(rep.keys()))
                            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))
                        else:
                            rep = {'ENV': '{:5.5}'.format(env_name), 'PW_ION': 'PW', 'LOW': '1.0', 'HIGH': '2.0'}
                            rep = dict((re.escape(k), v) for k, v in rep.items())
                            pattern = re.compile('|'.join(rep.keys()))
                            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))


"""
 Copies template MD parameter files to the simulation directories. Removes the 'WATER-', 'SOLVENT-' and charge flags. 
 Calls change_mdp().
"""


def copy_mdp(mdp_files, new_dir, env, ions):

    for mdp in mdp_files:

        new_path = Path(new_dir, mdp.name)
        try:
            shutil.copy(mdp, new_dir)
        except IOError as e:
            print('Copy operation failed: %s' % e)
        try:
            new_path.rename(Path(new_path.parent, new_path.stem.split('-')[1] + new_path.suffix))
            new_path = Path(new_dir, new_path.stem.split('-')[1] + new_path.suffix)
        except OSError:
            print('MDP-file renaming failed')

        try:
            change_mdp(new_path, env, ions)
        except OSError:
            print('Could not insert environment')


"""
 Handles the copying of the appropriate MD parameter template files to the simulation directories.
"""


def copy_files(env_path, mol_path, mdp_path, cluster_path, **kwargs):

    charged = re.compile(r'^CHARGED-')
    neutral = re.compile(r'^NEUTRAL-')

    for env in sorted(env_path):
        print(f'\nPreparing environment: {env.stem}\n')

        for idx, mol in enumerate(sorted(Path(mol_path).iterdir())):

            mol_name = mol.parts[-1]
            print(f'Copying molecule: {mol_name}')

            new_dir = Path(env, mol_name)

            copy_mols(mol, new_dir, mol_name)

            if 'contains_ions' not in kwargs:
                env_gro = env.stem.lower() + '.gro'
                gro_file = new_dir / env_gro
                ions = gro_file
            else:
                ions = kwargs['contains_ions'][str(env.parts[-1])][idx]

            if 'mol_charges' not in kwargs:
                mol_itp = mol_name + '.itp'
                itp_file = new_dir / mol_itp
                pos, neg = count_ions(itp_file)
                charge = pos - neg
            else:
                charge = kwargs['mol_charges'][mol_name]

            if env.stem == 'OCTANE':
                mdp_files = [path for path in Path(mdp_path).iterdir() if 'OCTANE-' in path.name]

            else:

                if charge == 0:
                    if env.stem == 'WATER':
                        mdp_files = [path for path in Path(mdp_path).iterdir() if 'WATER_NEUTRAL-' in path.name]
                    else:
                        mdp_files = [path for path in Path(mdp_path).iterdir() if re.search(neutral, path.name)]
                else:
                    if env.stem == 'WATER':
                        mdp_files = [path for path in Path(mdp_path).iterdir() if 'WATER_CHARGED-' in path.name]
                    else:
                        mdp_files = [path for path in Path(mdp_path).iterdir() if re.search(charged, path.name)]

            copy_mdp(mdp_files, new_dir, env, ions)
            copy_controls(cluster_path, new_dir, charge, mol_name, env.stem)


"""
 Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Automatic free energies of arbitrary small molecules.')
    parser.add_argument('-env', '--environment', type=Path, nargs='+', required=True,
                        help='Need: paths to one or more environment directories.')
    parser.add_argument('-mdp', '--parameters', type=Path, required=True,
                        help='Location of mdp options files. Prepend \"VACUUM-\" for small molecule energy '
                             'minimization and \"WATER-\" for simulation in water box to file names!')
    parser.add_argument('-mol', '--molecules', type=Path, required=True,
                        help='Need: location of molecule directories (parent directory).')
    parser.add_argument('-c', '--cluster', type=Path, required=True,
                        help='Need: location of setup files for running simulations on cluster.')

    args = parser.parse_args()

    copy_files(args.environment, args.molecules, args.parameters, args.cluster)
