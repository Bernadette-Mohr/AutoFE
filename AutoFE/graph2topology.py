import networkx as nx
import regex as re
import itertools
from pathlib import Path


"""
Reads forcefield file, finds [ atomtypes ] section, returns matches as dictionary (keys: atom/bead name, values: list of attributes)
"""


def get_forcefield(forcefield):
    ff_dict = {}

    """
    Use module regex: variable length lookbehind (?<=(; match.*\n)).
    Include '\n' to get multi-line matches (DOTALL), set '^', '$' to start and end of line (MULTILINE).
    Make non-greedy to return only first match: .*? (?)
    """
    rgx = re.compile(r'(?<=\[\s*atomtypes\s*\]).*?(?=\[\s*nonbond_params\s*\])', re.DOTALL | re.MULTILINE)

    with open(forcefield, 'rt') as ff_file:
        ff = ff_file.read()
        # Split the match by lines to get list of strings
        match = re.search(rgx, ff).group(0).split('\n')
        # Remove all blank lines and comments in the match
        atomtypes = [atom for atom in list(filter(None, match)) if not atom.startswith(';')]
        for s in atomtypes:
            ff_dict[s.split()[0]] = s.split()[1:]

    return ff_dict


"""
Prepend "S" to node name if node is in cycle node list, 
returns updated nodes dictionary.
"""


def set_cycles(nodes: object, cycles: object) -> object:
    for cn in cycles:
        nodes[cn] = 'S' + nodes[cn]

    return nodes


"""
Networkx utilities used to return a dictionary of node: bead/atom name, 
a unique list of all nodes that are part of a cycle
and a list of edges as (node, node) tuples
"""


def get_graph_attributes(graph):
    # Remove self looping edges
    graph.remove_edges_from(nx.selfloop_edges(graph))

    # Remove nodes that are not connected to the graph (degree = 0).
    # graph.remove_nodes_from(list(nx.isolates(graph)))
    # mapping = dict(zip(graph, range(0, graph.number_of_nodes() + 1)))
    # graph = nx.relabel_nodes(graph, mapping)

    # Get node labels from cleaned up graph.
    nodes = nx.get_node_attributes(graph, 'name')

    """
    All cycles in a graph can be constructed from the cycle basis - 
    all nodes that are part of any cycle in the graph are included in the cycle basis.
    """
    cycles = nx.cycle_basis(graph)
    cycle_nodes = []
    if not len(cycles) == 0:
        unique = set(itertools.chain.from_iterable(cycles))
        cycle_nodes.extend(list(unique))

    S_nodes = set_cycles(nodes, cycle_nodes)

    edges = list(graph.edges)

    degrees = {}
    for node, label in S_nodes.items():
        degrees[node] = graph.degree[node]

    return S_nodes, edges, degrees


"""
Charge state of beads are stored as sign in nodes dictionary values. 
Returns the corresponding float.
"""


def get_charge(name):
    charge = 0.0
    if '-' in name:
        charge = -1.0
    elif '+' in name:
        charge = 1.0

    return charge


"""
Lookup of bead/atom name in atomtypes dictionary, returns the mass 
of the bead/atom.
"""


def get_mass(bead, ff_dict):
    return ff_dict[bead][0]


"""
Writes the [atomtypes] section of the topology file.
"""


def write_atoms(nodes, itp_file, ff_dict):
    itp_file.write('[atoms]\n')
    itp_file.write('; i\ttype\tresnr\tresidue\tatom\tcgnr\tcharge\tmass\ttypeB\tchargeB\tmassB\n')

    beads = []
    charges = []

    for key, val in nodes.items():
        i = str(key + 1)
        bead = re.sub('\-$|\+$', '', val)
        beads.append(bead)
        charge = get_charge(val)
        charges.append(charge)
        itp_file.write(i + '\t' + bead + '\t' + str(1) + '\t' + 'MOL' + '\t' + bead + '\t' + i + '\t' + str(
            charge) + '\t' + get_mass(bead, ff_dict) + ';\n')
    itp_file.write('\n')

    return beads, charges


"""
Bonds are defined by edges in the graph, assigned different length and force constant depending on whether ring 
atoms/beads are involved.
"""


def write_bonds(edges, nodes, itp_file):
    itp_file.write('[bonds]\n')
    itp_file.write('; i\tj\tfunct\tlength\tforce.c.\n')

    bonds = []

    for edge in edges:

        i = str(edge[0] + 1)
        j = str(edge[1] + 1)
        # Normalize to boolean to cover XOR case
        i_ring = nodes[edge[0]].startswith('S')
        j_ring = nodes[edge[1]].startswith('S')

        if not i_ring and not j_ring:
            itp_file.write(i + '\t' + j + '\t' + '1' + '\t' + '0.37' + '\t' + '1250' + '\n')
            bonds.append((i, j, 0.37))
        elif i_ring ^ j_ring:
            itp_file.write(i + '\t' + j + '\t' + '1' + '\t' + '0.32' + '\t' + '1250' + '\n')
            bonds.append((i, j, 0.32))
        else:
            itp_file.write(i + '\t' + j + '\t' + '1' + '\t' + '0.25' + '\t' + '20000' + '\n')
            bonds.append((i, j, 0.25))

    itp_file.write('\n')

    return bonds


"""
Virtual sites needed for eg. for multiple aromatic rings
"""


def write_virtuals(itp_file):
    itp_file.write('[ virtual_sites3 ]\n')
    itp_file.write('; In-plane virtual sites defined by the three beads on the central ring\n')
    itp_file.write('; i\tj\tk\tl\tfunct\tjk\tkl\n')

    # TODO: add functionality if molecules can contain multiple rings

    itp_file.write('\n')


"""
Angles are set for every three consecutive atoms/beads in a chain 
up to the first connection to a ring.
"""


def write_angles(nodes, edges, degrees, itp_file):
    itp_file.write('[angles]\n')
    itp_file.write('; i\tj\tk\tfunct\tangle\tforce.c.\n')

    # TODO: extend for cases with longer chains
    start_nodes = [node for node, degree in degrees.items() if degree == 1]

    for node in start_nodes:
        # second_node can only contain one element
        second_node = [edge[1 - edge.index(node)] for edge in edges if node in edge][0]

        if degrees[second_node] == 2 and not nodes[second_node].startswith('S'):
            # will contain at least two elements
            next_nodes = [edge[1 - edge.index(second_node)] for edge in edges if second_node in edge]
            for third_node in next_nodes:
                if not third_node == node:
                    itp_file.write(str(node + 1) + '\t' + str(second_node + 1) + '\t' + str(
                        third_node + 1) + '\t' + '2' + '\t' + '180.0' + '\t' + '25.0\n')
        if degrees[second_node] == 3 and not nodes[second_node].startswith('S'):
            # will contain at least two elements
            next_nodes = [edge[1 - edge.index(second_node)] for edge in edges if second_node in edge]
            for third_node in next_nodes:
                if not third_node == node and not degrees[third_node] == 1:
                    itp_file.write(str(node + 1) + '\t' + str(second_node + 1) + '\t' + str(
                        third_node + 1) + '\t' + '2' + '\t' + '120.0' + '\t' + '25.0\n')
    itp_file.write('\n')


"""
Short-range interactions are excluded up to a distance of 2 in the graph 
and only for ring atoms/bonds.
"""


def get_exclusion_list(graph, nodes):
    nx.set_node_attributes(graph, nodes, 'name')

    excludes = {node: label for node, label in nodes.items() if not label.startswith('S')}
    graph.remove_nodes_from(excludes)
    rings = nx.get_node_attributes(graph, 'name')

    exclusions_list = []
    for node, label in rings.items():
        ego_dict = nx.get_node_attributes(nx.ego_graph(graph, node, radius=2), 'name')
        ex = [key + 1 for key in ego_dict.keys()]
        if len(ex) > 1:
            exclusions_list.append(ex)
        graph.remove_node(node)

    return exclusions_list


"""
Writes the exclusion list to the topology file.
"""


def write_exclusions(graph, nodes, itp_file):
    itp_file.write('[ exclusions ]\n')
    itp_file.write('; i\tj\tk\t...\n')

    exclusions_list = get_exclusion_list(graph, nodes)

    itp_file.write('\n'.join('\t'.join(map(str, row)) for row in exclusions_list))

    itp_file.write('\n')


"""
Handles creation of a simulation directory and the topology file for every molecule in the database.
"""


def create_topology(graphs_archive: object, forcefield: object) -> object:
    ff_dict = get_forcefield(forcefield)
    graphs = nx.read_gpickle(graphs_archive)
    root = ''

    """
    Number of digits of the number of graphs for formatting directory/file names. Replace with check if math.log10 
    returns an integer?
    """
    # n_digits = 0
    if len(graphs) == 100 or len(graphs) == 1000:
        n_digits = len(str(len(graphs))) - 1
    else:
        n_digits = len(str(len(graphs)))

    molecules_info = {}

    for idx, g in enumerate(graphs):

        mol = 'molecule_' + str(idx).zfill(n_digits)
        print(f'Processing {mol}')
        try:
            root = Path('molecules', mol)
            root.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            print(f'Topology for \"{mol}\" already exists, check path!')

        itp = mol + '.itp'
        itp_path = root / itp

        with open(itp_path, 'a+') as itp_file:

            itp_file.write(';;;; molecule ' + str(idx) + '\n')
            itp_file.write('\n')
            itp_file.write('[moleculetype]\n')
            itp_file.write('; molname\tnrexcl\n')
            itp_file.write('MOL\t4\n')
            itp_file.write('\n')

            G = graphs[idx]
            nodes, edges, degrees = get_graph_attributes(G)
            beads, charges = write_atoms(nodes, itp_file, ff_dict)
            bonds = write_bonds(edges, nodes, itp_file)
            write_virtuals(itp_file)
            write_angles(nodes, edges, degrees, itp_file)
            write_exclusions(G, nodes, itp_file)

        molecules_info[idx] = [beads, charges, bonds]

        top = mol + '.top'
        top_path = root / top

        with open(top_path, 'a+') as top_file:

            top_file.write('#include \"' + str(forcefield.absolute()) + '\"\n')
            top_file.write(f'#include \"{mol}.itp\"\n')  # {mol}/
            top_file.write('\n')
            top_file.write('[ system ]\n')
            top_file.write(f'{mol} in vacuum\n')
            top_file.write('\n')
            top_file.write('[ molecules ]\n')
            top_file.write('; name\tnumber\n')
            top_file.write('MOL\t1\n')

    return root.parent, molecules_info


"""
Added functionality for running module standalone
"""
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create Gromacs topologies from molecule graphs.')
    parser.add_argument('-mol', '--molecules', type=Path, required=True,
                        help='Need: location of molecule graphs archive')
    parser.add_argument('-ff', '--forcefield', type=Path, required=True, help='Location of main forcefield file.')

    args = parser.parse_args()

    create_topology(graphs_archive=args.molecules, forcefield=args.forcefield)
