import os, sys
import gromacs
from pathlib import Path
import fileinput
from contextlib import closing
import regex as re


# TODO: find way to handle context manager without needing os. Maybe contextlib contextmanager?
class cd:
    """Context manager for safely changing the current working directory"""

    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


"""
 For minimization round with increasing timesteps. Copies the template to the molecule directory, replaces the placeholder variable with the required time step.
"""


def set_timestep(mdp_file, dt):
    template = open(mdp_file, 'rt').read()
    param_file = Path(f'VACUUM-3.{dt}.mdp')

    with open(param_file, 'wt') as output:
        output.write(template)

    with closing(fileinput.FileInput(param_file, inplace=True)) as pf:
        for line in pf:
            sys.stdout.write(line.replace('DELTA_T', str(dt)))

    return param_file


"""
 If Fmax of the molecule doesn't get minimized to below 1 kcal/mol within 50 attempts, the threshold Ftol is increased to 10 kcal/mol.
"""


def adjust_Ftol(mdp_file, mol=None):
    template = open(mdp_file, 'rt').read()
    new_file = ''

    if mol is not None:
        new_file = mol / mdp_file.name
    else:
        new_file = Path(mdp_file.name)

    print('Adjusted mdp', new_file)

    with open(new_file, 'wt') as output:
        output.write(template)

    with closing(fileinput.FileInput(new_file, inplace=True)) as pf:
        for line in pf:
            sys.stdout.write(line.replace('emtol                    = 1', 'emtol                    = 10'))

    return new_file


"""
 Check the Gromacs log file for messages indicating that the minimization was not successful. Boolean.
"""


def check_if_minimized(logfile):
    rgx = re.compile(
        r'(Steepest Descents did not converge to Fmax <)|(Steepest Descents converged to machine precision)',
        re.MULTILINE)

    with open(logfile, 'rt') as log_file:

        log = log_file.read()
        found = re.search(rgx, log)

        if found is not None:
            return True
        else:
            return False


"""
 Handles the minimization in vacuum/gas phase of the molecules with randomly placed particles.
"""


def run_minimization(mols_path, mdp_path):
    """
     Using absolute paths here allows for later changing into the molecule subdirectories for running Gromacs.
    """
    mdp_files = [path.absolute() for path in sorted(Path(mdp_path).iterdir()) if 'VACUUM-' in path.name]
    timesteps = [0.0002, 0.0005, 0.001, 0.002]

    for mol in sorted(mols_path.iterdir()):

        name = mol.stem

        with cd(mol):

            """
             Attribute 'gro_file' will be overwritten with the current output of each minimization round
            """
            gro_file = Path(name + '.gro')
            print('Starting', gro_file)
            for mdp_file in mdp_files:

                number = mdp_file.name.split('-')[1].split('.')[0]
                top_file = Path(name + '.top')
                tpr_file = Path(name + '-' + number + '.tpr')

                if mdp_file.name.endswith('-3.mdp'):

                    for dt in timesteps:

                        mdp = set_timestep(mdp_file, dt)

                        gromacs.grompp(f=mdp, c=gro_file, p=top_file, o=tpr_file)
                        gromacs.mdrun(v=False, deffnm=name + '-' + number)

                        gro_file = Path(name + '-' + number + '.gro')

                        """
                         counter to emulate do{...}while() construct for relaxing Ftol after 50 failed minimization 
                         attempts
                        """
                        if_counter = 0

                        while check_if_minimized(name + '-' + number + '.log'):

                            if_counter += 1

                            if if_counter == 50:
                                mdp = adjust_Ftol(mdp)

                            gromacs.grompp(f=mdp, c=gro_file, p=top_file, o=tpr_file)
                            gromacs.mdrun(v=False, deffnm=name + '-' + number)

                            gro_file = Path(name + '-' + number + '.gro')

                else:

                    gromacs.grompp(f=mdp_file, c=gro_file, p=top_file, o=tpr_file)
                    gromacs.mdrun(v=False, deffnm=name + '-' + number)

                    gro_file = Path(name + '-' + number + '.gro')

                    else_counter = 0

                    while check_if_minimized(name + '-' + number + '.log'):

                        """
                         counter to emulate do{...}while() construct for relaxing Ftol after 50 failed minimization 
                         attempts
                        """
                        else_counter += 1

                        if else_counter == 50:
                            mdp_file = adjust_Ftol(mdp_file, mol)

                        gromacs.grompp(f=mdp_file, c=gro_file, p=top_file, o=tpr_file)
                        gromacs.mdrun(v=False, deffnm=name + '-' + number)

                        gro_file = Path(name + '-' + number + '.gro')

        # All temporary run and Gromacs backup files are deleted.
        for backup in mol.glob('#*'):
            backup.unlink()
        for tmp_mdp in mol.glob('*.mdp'):
            tmp_mdp.unlink()
        # TODO: provide option to automatically delete all output files from intermediate minimization steps
        #  (keep only input, step 4)


"""
Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Automatic free energies of arbitrary small molecules.')
    parser.add_argument('-mol', '--molecules', type=Path, required=True,
                        help='Need: location of molecule directories (parent directory).')
    parser.add_argument('-mdp', '--parameters', type=Path, required=True,
                        help='Location of mdp options files. Prepend \"VACUUM-\" for small molecule energy '
                             'minimization and \"WATER-\" for simulation in water box.')

    args = parser.parse_args()

    run_minimization(args.molecules, args.parameters)
