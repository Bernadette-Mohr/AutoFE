import sys
import regex as re
import random
import shutil
import fileinput
from contextlib import closing
from pathlib import Path

"""
 ONLY CALLED IF MODULE IS EXECUTED STANDALONE
"""


def count_ions(mol_itp):
    #
    rgx = re.compile(r'(?<=\[atoms\]\n).*?(?=^$)', re.DOTALL | re.MULTILINE)

    with open(mol_itp, 'rt') as itp_file:
        #
        itp = itp_file.read()
        match = re.search(rgx, itp).group(0).split('\n')
        atoms = [atom for atom in list(filter(None, match)) if not atom.startswith(';')]
        charges = [charge.split('\t')[6] for charge in atoms]
        negatives = [n for n in charges if re.search(r'^\-1.0', n)]
        positives = [p for p in charges if re.search(r'(^1.0)|(\+1.0)', p)]

        return len(positives), len(negatives)


"""
 Get the number of excess charges (not internally equalized) in the molecule.
 Returns number of excess charges and flag if it is positive or negative.
"""


def get_charge_difference(n_positives, n_negatives):
    #
    if n_positives >= n_negatives:
        return n_positives - n_negatives, 'positive'
    else:
        return n_negatives - n_positives, 'negative'


"""
 Copy environment files into environment/molecule_* subdirectories.
"""


def make_envdir(new_dir, env_file):
    #
    files = [env_file, env_file.with_suffix('.top')]

    for f in files:
        try:
            shutil.copy(f, new_dir)
        except IOError as e:
            print('Copy operation failed: %s' % e)


"""
 Remove the required number from the ion count in the gromacs 
 *.top file or add the inserted number of counterions.
"""


def adapt_topology(n_charges, ion, env_gro, new_dir, no_ion=False):
    #
    top = env_gro.with_suffix('.top').name
    new_top = Path(new_dir) / top
    rgx = re.compile(ion)
    found = False
    with closing(fileinput.FileInput(new_top, inplace=True)) as new_file:
        file_content = fileinput.input(new_top)
        lines = [line for line in file_content]
        last = lines[-1]
        for line in new_file:
            if not no_ion:
                if re.search(rgx, line):
                    sys.stdout.write(line.split()[0] + '\t' + str(int(line.split()[1]) - n_charges) + '\n')
                    found = True
                elif 'PW' in line:
                    sys.stdout.write(line.split()[0] + '\t' + str(int(line.split()[1]) + n_charges) + '\n')
                else:
                    sys.stdout.write(line)
            else:
                if re.search(rgx, line):
                    sys.stdout.write(line.split()[0] + '\t' + str(int(line.split()[1]) + n_charges) + '\n')
                    found = True
                elif 'PW' in line:
                    sys.stdout.write(line.split()[0] + '\t' + str(int(line.split()[1]) - n_charges) + '\n')
                else:
                    sys.stdout.write(line)
            if line == last and not found:
                print(ion + '\t' + str(n_charges))


"""
 Find all occurrences of the required ion (as regex pattern rgx) in the environment gro file. Typecast to set 
 automatically removes all duplicates (3-body ion: 3 times ion index). Returns none if ion is not present in gro file.
"""


def contains_mols(gro, rgx=None):
    if rgx is not None:
        return set(re.findall(rgx, gro))
    else:
        print('No ions in environment!')


"""
 Picks a number n of ions randomly out of all ions of the required kind in the gro file.
"""


def pick_random(match, n):
    return random.sample(match, n)


"""
 Rename ions to water, water to ions.
 Regex pattern using regex syntax cannot be used as dictionary keys, only plain strings. 
 So renaming water to ion catches all cases for the first pattern, results have to be corrected.
"""


def rename_molecules(mols, line, rename_dict):
    #
    pattern = re.compile('|'.join(rename_dict.keys()))
    for mol in mols:
        if mol in line:
            line = pattern.sub(lambda m: rename_dict[m.group()], line)
            if 'NACP' in line:
                line = line.replace('NACP', ' NAP')
            if 'NACM' in line:
                line = line.replace('NACM', ' NAM')
            if 'CLCP' in line:
                line = line.replace('CLCP', ' CLP')
            if 'CLCM' in line:
                line = line.replace('CLCM', ' CLM')
    return line


"""
 Checks if ions with the required charge are present in the environment. If yes, it renames the required number to 
 water. If no, it adds the required number of counterions. Returns True if there are any ions in the resulting 
 coordinate file, False otherwise.
"""


def adjust_charges(n_charges, ion, env_gro, new_dir):

    new_gro = new_dir / env_gro.name

    pos = re.compile(r'\d+PNA')
    neg = re.compile(r'\d+PCL')
    wat = re.compile(r'\d+PW')

    water = ['PW ', '  W', ' WP', ' WM']
    sodium = ['PNA', 'NAC', 'NAP', 'NAM']
    chlorine = ['PCL', 'CLC', 'CLP', 'CLM']

    if n_charges == 0:
        return new_gro, pos, neg
    else:
        print(f'Need to adjust for {n_charges} {ion} charges')

    random_ions, random_solvent = [], []

    with open(env_gro, 'rt') as gro_file:
        #
        gro = gro_file.read()

        if 'PW' not in gro:
            return new_gro, pos, neg
        else:
            if ion == 'positive':
                match = contains_mols(gro, pos)
                ion_type = 'PNA'
                ion_dict = dict(zip(sodium, water))
                water_dict = dict(zip(water, chlorine))
            else:
                match = contains_mols(gro, neg)
                ion_type = 'PCL'
                ion_dict = dict(zip(chlorine, water))
                water_dict = dict(zip(water, sodium))

            if match:
                out_ion = ion.capitalize()
                print(f'{out_ion} ions in environment')
                random_ions = pick_random(match, n_charges)
                adapt_topology(n_charges, ion_type, env_gro, new_dir)
            else:
                print(f'NO {ion} ions in environment')
                solvent = contains_mols(gro, wat)
                random_solvent = pick_random(solvent, n_charges)

                no_ion = True
                if ion_type == 'PNA':
                    new_ion = 'PCL'
                    adapt_topology(n_charges, new_ion, env_gro, new_dir, no_ion)
                else:
                    new_ion = 'PNA'
                    adapt_topology(n_charges, new_ion, env_gro, new_dir, no_ion)

            with closing(fileinput.FileInput(new_gro, inplace=True)) as new_file:
                for line in new_file:
                    if random_ions:
                        sys.stdout.write(rename_molecules(random_ions, line, ion_dict))
                    else:
                        sys.stdout.write(rename_molecules(random_solvent, line, water_dict))

    return new_gro, pos, neg


"""
 Checks if there is any kind of ion in the adapted environment.gro file for adapting the coupling groups in the mdp 
 files. Boolean
"""


def has_ions(new_gro, pos, neg):
    #
    ions = False

    with open(new_gro, 'rt') as new_file:
        new_content = new_file.read()
        if contains_mols(new_content, pos) or contains_mols(new_content, neg):
            ions = True

    return ions


"""
 Works only with polarizable water and ion forcefield with monovalent ions.
 Works for apolar solvents.
 Takes coordinate file of environment and topology of small molecule as arguments. 
 Gets the number of positive and negative ions present in the molecule. 
 Checks environment for presence of ions. If yes, required number of respective ion is renamed to water. If not, 
 required number of water molecules is renamed to respective ion.
"""


def equilibrate(environments, molecules, molecules_info=None):
    
    root = molecules.parent

    env_paths = []
    contains_ions = dict()
    mol_charge = dict()
    
    for env_gro in sorted(environments.glob('*.gro')):

        env = env_gro.stem.upper()
        new_path = Path(root, env)
        env_paths.append(new_path)
        print(f'\nEnvironment: {env}\n')

        contains_ions[env] = dict()

        for idx, mol_path in enumerate(sorted(molecules.iterdir())):

            mol = mol_path.parts[-1]
            itp_name = mol + '.itp'
            mol_itp = mol_path / itp_name
            try:
                new_dir = Path(new_path, mol)
                new_dir.mkdir(parents=True, exist_ok=False)
            except FileExistsError:
                print(f'Envionments for \"{mol}\" already equilibrated, check path!')

            print(f'\nProcessing {mol}\n')

            if molecules_info is None:
                n_positives, n_negatives = count_ions(mol_itp)
            else:
                n_positives = len([p for p in molecules_info[idx][1] if p > 0])
                n_negatives = len([n for n in molecules_info[idx][1] if n < 0])

            n_charges, ion = get_charge_difference(n_positives, n_negatives)

            if n_charges == 0:
                print('Molecule charge at equilibrium!')
                make_envdir(new_dir, env_gro)
                new_gro, pos, neg = adjust_charges(n_charges, ion, env_gro, new_dir)
                mol_charge[mol] = n_charges
                contains_ions[env][idx] = has_ions(new_gro, pos, neg)
            else:
                make_envdir(new_dir, env_gro)
                new_gro, pos, neg = adjust_charges(n_charges, ion, env_gro, new_dir)
                mol_charge[mol] = n_charges
                contains_ions[env][idx] = has_ions(new_gro, pos, neg)

    return env_paths, mol_charge, contains_ions


"""
Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Automatic free energies of arbitrary small molecules.')
    parser.add_argument('-env', '--environment', type=Path, required=True, help='Need: location of the Gromacs input '
                                                                                'files for the simulation environments.')
    parser.add_argument('-mol', '--molecules', type=Path, required=True, help='Need: location of molecule directories '
                                                                              '(parent directory).')

    args = parser.parse_args()

    equilibrate(args.environment, args.molecules)
