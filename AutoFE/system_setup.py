import sys
import gromacs as gmx
from pathlib import Path
import regex as re
from contextlib import closing
import fileinput


"""
 Create the initial coordinate file of the whole system by placing the small molecule at the appropriate position in 
 the environment.
"""


def make_init(env, child, mol):
    print('Creating initial system coordinates')

    env_name = env.stem.lower() + '.gro'
    mol_name = mol + '.gro'
    mol_tmp = mol + '-2.gro'
    env_pdb_tmp = env.stem.lower() + '.pdb'
    mol_pdb_tmp = mol + '.pdb'

    init = child / 'init.gro'

    env_gro = child / env_name
    mol_gro = child / mol_name
    mol_center = child / mol_tmp
    env_pdb = child / env_pdb_tmp
    mol_pdb = child / mol_pdb_tmp

    size = env_gro.read_text().rstrip('\n').split('\n')[-1].split()

    """
     In the water or octane box, the molecule is placed at the box center.
    """
    if env.stem == 'WATER' or env.stem == 'OCTANE':
        center = [float(x) / 2 for x in size]
        gmx.editconf(f=mol_gro, o=mol_center, center=center)
        gmx.editconf(f=mol_center, o=mol_pdb)
        gmx.editconf(f=env_gro, o=env_pdb)
    else:
        """
         With a lipid bilayer, the molecule is placed at the interface region of the membrane.
        """
        center = [float(size[0]) / 2, float(size[1]) / 2, (float(size[2]) / 2) - 2]
        gmx.editconf(f=mol_gro, o=mol_center, center=center)
        gmx.editconf(f=mol_center, o=mol_pdb)
        gmx.editconf(f=env_gro, o=env_pdb)

    with mol_pdb.open('a') as tmp:
        tmp.write(env_pdb.read_text())

    with closing(fileinput.FileInput(mol_pdb, inplace=True)) as tmpfile:
        for line in tmpfile:
            if 'ATOM' in line:
                sys.stdout.write(line)

    gmx.editconf(f=mol_pdb, o=init, box=size)

    with closing(fileinput.FileInput(init, inplace=True)) as init_file:
        for line in init_file:
            rep = {'PNA': 'ION', 'PCL': 'ION'}
            rep = dict((re.escape(k), v) for k, v in rep.items())
            pattern = re.compile('|'.join(rep.keys()))
            sys.stdout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))

    for ext in ('#*', '*.pdb', '*-2.gro'):
        for _file in child.glob(ext):
            _file.unlink()


"""
 Add the molecule topology file and the number of molecules present to the environment.top file. Rename the resulting 
 file to 'system.top'.
"""


def make_topology(env, child, mol):
    print(f'Adding {mol} to top file')

    env_name = env.stem.lower() + '.top'
    env_top = child / env_name
    system = child / 'system.top'

    rgx = re.compile(r'(?i)include.*\s')

    with closing(fileinput.FileInput(env_top, inplace=True)) as top:
        #

        match = False
        skip = False

        for line in top:

            if re.search(rgx, line) and not match:
                match = True

            if line.isspace() and match:
                print(line.replace(line, '; Atom types of specific molecule\n#include \"' + mol + '.itp\"', ))
                match = False

            if not skip:
                if '[ system ]' in line:
                    print(line.replace(line, line + mol + '-' + env.stem + '\n'))
                    skip = True
                else:
                    if len(line.split()) == 2 and not re.search(rgx, line):
                        sys.stdout.write(line.split(None, 1)[0] + '\t' + line.split(None, 1)[1])
                    else:
                        sys.stdout.write(line)

            if skip and line.isspace():
                skip = False

            if 'number' in line:
                print('MOL\t1')

    env_top.rename(system)


"""
 Generating the Gromacs system index file with the coupling groups used in the simulation.
"""


def make_index(env, child, mol):
    print(f'Creating index file for system {env}-{mol}')

    init = child / 'init.gro'
    idx = child / 'index.ndx'

    content = init.read_text()

    if env.stem == 'WATER' or env.stem == 'OCTANE':
        if 'ION' in content:
            gmx.make_ndx(f=init, o=idx, input=('3 | 4', 'q'))
        else:
            gmx.make_ndx(f=init, o=idx, input=('q'))
    else:
        if 'ION' in content:
            gmx.make_ndx(f=init, o=idx, input=('2 | 3', '4 | 5', 'q'))
        else:
            gmx.make_ndx(f=init, o=idx, input=('2 | 3', 'q'))

    env_file = env.stem.lower() + '.gro'
    for ext in ('#*', 'molecule_*.gro', 'molecule_*.top', env_file):
        for _file in child.glob(ext):
            _file.unlink()


"""
 For every subdirectory in the environment directories, the initial systems are created from the environment and 
 molecule files.
"""


def setup_system(environments):
    #
    for env in sorted(environments):
        for child in sorted(env.glob('molecule_*')):
            #
            mol_name = str(child.stem)

            make_init(env, child, mol_name)
            make_topology(env, child, mol_name)
            make_index(env, child, mol_name)


"""
 Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='System setup for simulating arbitrary small molecules in different environments.')
    parser.add_argument('-env', '--environment', type=Path, nargs='+', required=True,
                        help='Need: paths to one or more environment directories.')

    args = parser.parse_args()

    setup_system(args.environment)
