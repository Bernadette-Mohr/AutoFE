import sys
from pathlib import Path
import shutil
from contextlib import closing
import fileinput
from equilibrate_charge import count_ions

"""
 Create names for the lambda subdirectories containing environment, molecule and lambda step for future reference.
"""


def set_dirname(mol_path, i, n_digits):
    #
    return mol_path.parts[-2] + '-' + str(mol_path.stem) + '-' + str(i).zfill(n_digits)


"""
 Create the subdirectory for each lambda step in the molecule simulation directories.
 To be removed in the future since pathlib.Path.mkdirs() creates all missing parts of a path.
"""


def make_dir(mol_path, dirname):
    #
    path = mol_path / dirname
    path.mkdir(parents=False, exist_ok=False)
    return path


"""
 Set the lambda step in the MD parameter files by replacing the placeholder 'LAMBDA' in the template with the 
 corresponding integer.
"""


def set_lambda(filename, i):
    #
    with closing(fileinput.FileInput(filename, inplace=True)) as mdp_file:
        for line in mdp_file:
            sys.stdout.write(line.replace('LAMBDA', str(i)))


"""
 Copy all system and MD parameter files to the lambda subdirectories and adjust their names to the corresponding 
 lambda step.
"""


def copy(mol_path, new_path, i, n_digits):
    #
    if i == 0:
        init = mol_path / 'init.gro'
        first_init = init.stem + '-' + str(i).zfill(n_digits) + init.suffix
        new_init = mol_path / first_init
        init.rename(new_init)
        shutil.copy(new_init, new_path)

    local_files = sorted([path for path in mol_path.glob('*') if path.suffix and not path.suffix == '.gro'])
    for lf in local_files:
        if lf.is_file():
            shutil.copy(lf, new_path)

    for mdp in new_path.glob('*.mdp'):
        new_name = mdp.stem + '-' + str(i).zfill(n_digits) + mdp.suffix
        new_file = new_path / new_name
        mdp.rename(new_file)
        set_lambda(new_file, i)


"""
 For every molecule subdirectory in the environment simulation directories, create a sub-subdirectory for each lambda 
 step of the FE simulation.
"""


def create(environments, forcefield_file, charges=None):
    #
    forcefield = forcefield_file.parent

    for env in sorted(environments):
        print(f'\nPreparing lambda directories in {env.stem}\n')

        subdirs = sorted(env.glob('molecule_*'))

        if env.stem == 'OCTANE':
            topologies = [forcefield / 'martini_v2.2RP_5T.itp', forcefield / 'martini_v2.0_solvents.itp']
        elif env.stem == 'WATER':
            topologies = [forcefield / 'martini_v2.2RP_5T.itp', forcefield / 'ions.itp']
        else:
            topologies = [forcefield / 'martini_v2.2RP_5T.itp',
                          forcefield / 'martini_v2.0_lipids_all_201506.itp', forcefield / 'ions.itp']

        """
         Number of digits of the number of graphs for formatting directory/file names. Replace with check if math.log10 
         returns an integer?
        """
        length = len(subdirs)

        if length == 100 or length == 1000:
            n_digits = len(str(length)) - 1
        else:
            if length < 10:
                n_digits = 2
            else:
                n_digits = len(str(length))

        for mol in subdirs:
            print(f'... {mol.stem}')

            if charges is None:
                mol_itp = mol.stem + '.itp'
                itp_file = mol / mol_itp
                pos, neg = count_ions(itp_file)
                charge = pos - neg
            else:
                charge = charges[mol.stem]

            if charge != 0:
                if env.stem == 'OCTANE':
                    for i in range(0, 40):
                        new_path = make_dir(mol, set_dirname(mol, i, n_digits))
                        copy(mol, new_path, i, n_digits)
                else:
                    for i in range(0, 80):
                        new_path = make_dir(mol, set_dirname(mol, i, n_digits))
                        copy(mol, new_path, i, n_digits)
            else:
                for i in range(0, 40):
                    new_path = make_dir(mol, set_dirname(mol, i, n_digits))
                    copy(mol, new_path, i, n_digits)

            for _lambda in sorted(mol.iterdir()):
                if _lambda.is_dir():
                    for top in topologies:
                        shutil.copy(top, _lambda)


"""
 Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate subdirectories for all lambda steps in a free energy '
                                                 'simulation and copy the relevant files.')
    parser.add_argument('-env', '--environment', type=Path, nargs='+', required=True,
                        help='Need: paths to one or more environment directories.')
    parser.add_argument('-ff', '--forcefield', type=Path, required=True,
                        help='Location of main forcefield file to be used for molecule generation.')

    args = parser.parse_args()

    create(args.environment, args.forcefield)
