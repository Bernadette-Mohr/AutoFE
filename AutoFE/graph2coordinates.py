import numpy as np
import regex as re
from pathlib import Path

"""
ONLY CALLED IF MODULE IS EXECUTED STANDALONE
Get the beads from the topology file.
Returns list of the beads.
"""


def get_atom_dict(top_file):
    rgx = re.compile(r'(?<=\[atoms\].*\n).*?(?=^$)', re.DOTALL | re.MULTILINE)
    match = re.search(rgx, top_file).group(0).split('\n')
    atomtypes = [atom for atom in list(filter(None, match)) if not atom.startswith(';')]
    atoms = dict(zip([atom.split('\t')[0] for atom in atomtypes], [atom.split('\t')[1] for atom in atomtypes]))

    return atoms


"""
ONLY CALLED IF MODULE IS EXECUTED STANDALONE
Reads the topology file of the small molecule to get several attributes. Returns a list of included atoms, the sum of all included bond lengths, and the length of the shortest bond in the molecule + the bead diameter.
"""


def read_topology(top) -> object:
    with open(top, 'rt') as top_file:

        tf = top_file.read()
        atoms = get_atom_dict(tf)
        rgx = re.compile(r'(?<=\[bonds\]\n).*?(?=^$)', re.DOTALL | re.MULTILINE)
        match = re.search(rgx, tf).group(0).split('\n')
        bonds = [atom for atom in list(filter(None, match)) if not atom.startswith(';')]
        # Scaling bond lengths to make the initial sphere as small as possible
        if len(atoms) < 4:
            bonds_length = sum([float(bond.split('\t')[3]) for bond in bonds]) * 1.1
        else:
            bonds_length = sum([float(bond.split('\t')[3]) for bond in bonds]) * 0.85
        # add martini bead diameter (0.526 nm) to the minimum distance to prevent collisions.
        if len(bonds) is 0:
            shortest_bond = 0.526
        else:
            shortest_bond = min([float(bond.split('\t')[3]) for bond in bonds]) + 0.526

        return atoms, bonds_length, shortest_bond


"""
Creates a random unit vector in 3D for randomly placing a particle within a sphere of radius: sum of all bonds. Returns the x, y and z coordinates of that vector.
"""


def random_unit_vector(bonds_length):
    phi = np.random.uniform(0, np.pi * 2)
    costheta = np.random.uniform(-1, 1)

    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return x, y, z


"""
Picks a random factor between the length of the shortes bond and the sum of all bonds as a scaling factor for placing the new particle in 
the direction set by the random unit vector in the center of the box.
Retund the x, y and z coordinates of the new particle.
"""


def translate(bonds_length, shortest_bond):
    d = np.random.uniform(shortest_bond, bonds_length)

    x, y, z = random_unit_vector(bonds_length)

    u_x = (x + 2.5) * d
    u_y = (y + 2.5) * d
    u_z = (z + 2.5) * d

    return u_x, u_y, u_z


"""
Calculates the distance between two points in 3D space for collision check.
"""


def calculate_distance(x, y, z, test_x, test_y, test_z):
    distance = np.sqrt(np.square(test_x - x) + np.square(test_y - y) + np.square(test_z - z))

    return distance


"""
Checks if the new particle would be positioned with a distance of at least shortest_bond to every other particle, but still within the sphere radius of bonds_length. Boolean.
"""


def collision_check(positions, test_x, test_y, test_z, bonds_length, shortest_bond):
    lower = shortest_bond
    upper = bonds_length

    for key, val in positions.items():
        x = val[1][0]
        y = val[1][1]
        z = val[1][2]

        if lower <= calculate_distance(x, y, z, test_x, test_y, test_z) <= upper:
            return False
        else:
            return True


"""
Handles the positioning of all the beads in the molecule inside a spherical volume around the center of the simulation box. Returns the 
accepted positions for all beads as a dictionary with bead index as key and the Gromacs coordinate file line as value.
"""


def randomize_positions(atoms: object, bonds_length: object, shortest_bond: object) -> object:
    positions = {}

    for idx, atom in atoms.items():
        """
        The first bead is simply placed in the center of the simulation box.
        """
        if idx == '1':
            coords = (atom, [round(2.500, 3), round(2.500, 3), round(2.500, 3)])
            positions.update({idx: coords})
            print(f'{atom}')
        else:
            """
            All subsequent beads have to be placed inside a sphere with radius: sum of lengths of all bonds in the 
            molecule.
            """
            x, y, z = translate(bonds_length, shortest_bond)
            iteration = 0
            """
            As long as the next particle would be placed outside of the sphere or closer than the minimum distance to 
            any previously placed particle, the coordinates will be rejected and another placement attempt will be made.
            """
            while collision_check(positions, x, y, z, bonds_length, shortest_bond):
                iteration += 1

                x, y, z = translate(bonds_length, shortest_bond)

            print(f'{atom}, # iterations: ', iteration)

            coords = (atom, [round(x, 3), round(y, 3), round(z, 3)])
            positions.update({idx: coords})

    return positions


"""
 Handles generation of and output of the final accepted positions for all the beads in the molecule in Gromac coordinate file format.
"""


def write_coordinates(top, info=None):
    print(top.stem)

    """
    atoms: dict of all the particle indices (keys) and names (values) in the molecule
    bonds_length: sum of all bonds in the molecule to determine volume of sphere (float)
    shortest_bond: minimum distance all the particles have to keep to each other (float)
    positions: positions for all particles in the molecule meeting the imposed conditions.
    """
    # atoms = {}
    # bonds_length, shortest_bond = 0.0, 0.0

    if info is None:
        atoms, bonds_length, shortest_bond = read_topology(top)
    else:
        atoms = {str(k + 1): v for k, v in enumerate(info[0])}
        # bonds_length = sum([bond[2] for bond in info[2]])*0.85
        # shortest_bond = min([bond[2] for bond in info[2]]) + 0.526
        # Scaling bond lengths to make the initial sphere as small as possible
        if len(atoms) < 4:
            bonds_length = sum([bond[2] for bond in info[2]]) * 2.95
        else:
            bonds_length = sum([bond[2] for bond in info[2]]) * 0.95
        # add martini bead diameter (0.526 nm) to the minimum distance to help prevent collisions.
        if len(info[2]) is 0:
            shortest_bond = 0.526
        else:
            shortest_bond = min([bond[2] for bond in info[2]]) + 0.526

    positions = randomize_positions(atoms, bonds_length, shortest_bond)

    gro_name = top.stem + '.gro'
    path = top.parent / gro_name

    with open(path, 'a+') as gro_file:

        gro_file.write(top.stem + ': random-generated initial coordinates\n')
        gro_file.write(max(list(atoms.keys())).rjust(5) + '\n')

        for key, val in positions.items():
            x = str(val[1][0])
            y = str(val[1][1])
            z = str(val[1][2])

            gro_file.write(
                key.rjust(5) + 'MOL'.rjust(5) + val[0].rjust(5) + key.rjust(5) + x.rjust(8) + y.rjust(8) + z.rjust(8) + '\n')

        gro_file.write('   5.00000   5.00000   5.00000\n')


"""
 Iterates over all the processed molecules to create an initial coordinate file with arbitrary positioning of the 
 particles.
"""


def create_coordinates(topologies, molecules_info=None):
    for idx, top in enumerate(sorted(topologies.glob('**/*.itp'))):

        if top.is_file():

            if molecules_info is None:
                write_coordinates(top)
            else:
                write_coordinates(top, molecules_info[idx])


"""
Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate random coordinates for all beads in arbitrary '
                                                 'small molecules.')
    parser.add_argument('-mol', '--molecules', type=Path, required=True, help='Need: location of molecule directories '
                                                                              '(parent directory).')

    args = parser.parse_args()

    create_coordinates(args.molecules)
