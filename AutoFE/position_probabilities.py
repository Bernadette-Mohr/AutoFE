import argparse
from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import sys


def read_data(res_dirs, mol_name):
    rgx = re.compile(r'prod-\d{2,3}_pullx\.xvg')
    xvg_files = list()

    for res_dir in res_dirs:
        [xvg_files.append(_file) for _file in sorted(res_dir.iterdir()) if re.match(rgx, _file.name)]

    df = pd.DataFrame()
    t, x = [], []
    xvg = xvg_files[-1]

    with open(xvg, 'r') as infile:
        content = [line.strip() for line in infile.readlines() if not line.startswith('#') and not line.startswith('@')]
        for line in content:
            parts = line.split('\t')
            t.append(float(parts[0]))
            if parts[1] == parts[2]:
                x.append(float(parts[1]))
            else:
                print('Something wrong here!')
                break
    df['time'] = t
    df[mol_name] = x
    df = df.set_index('time')

    return df


def is_perfect_square(n_mols):

    root = np.sqrt(n_mols)

    if int(root + 0.5) ** 2 == n_mols:
        return [int(root)]
    else:
        return []


def perfect_squares(n_graphs):

    squares = list()

    for i in range(1, n_graphs + 1):
        #
        if np.sqrt(i) == np.trunc(np.sqrt(i)):
            squares.append(i)

    return np.asarray(squares, dtype=np.int16)


def set_diameter(n_mols, all_squares):

    sqrt_result = is_perfect_square(n_mols)

    if sqrt_result:
        n_rows = sqrt_result[0]
        n_cols = sqrt_result[0]
    else:
        n_rows = int(np.sqrt(all_squares[all_squares < n_mols].max()))
        n_cols = int(np.ceil(n_mols / n_rows))

    return n_rows, n_cols


def split_mols(mols, n):
    for i in range(0, len(mols), n):
        yield mols[i:i+n]


def set_figsize(n_rows, n_cols):

    if n_rows == n_cols:
        if n_rows == 1 and n_cols == 1:
            SIZE = (6, 5)
        elif n_cols < 4:
            SIZE = (15, 15)
        else:
            SIZE = (20, 20)
    elif n_cols < 5:
        SIZE = (15, 7)
    else:
        SIZE = (18, 10)

    return SIZE


def plot_molecule(mol, stats, midplane_cut, water_cut, X, ax):
    ax.axvspan(midplane_cut, water_cut, alpha=0.2, color='g')
    text_x = (water_cut - midplane_cut) * 0.5 + midplane_cut
    ax.text(text_x, 6, 'interfacial', horizontalalignment='center')

    ax.plot(X, stats[1], color='k', label=f'from midplane', linewidth=1.5, linestyle='-')
    ax.plot(X, stats[2], color='k', label=f'from water', linewidth=1.5, linestyle='--')

    ax.plot(stats[0][0], stats[0][1], 'x', color='r')
    ax.axvline(x=stats[0][0], ymax=0.5, color='r', linestyle='--')

    ax.set_xticks(X)
    ax.set_xlim(0.5, 10.5)
    ax.legend()
    ax.set_ylabel(r"$\sum\ \mathrm{probability} $")
    ax.set_xlabel('bins')


def plot_molecules(stats, env_path, midplane_cut, water_cut, X):

    n_mols = len(stats.keys())
    all_squares = perfect_squares(n_mols)

    if n_mols > 10:

        batches = split_mols(list(stats.keys()), 10)
        for batch in batches:
            n = len(batch)
            n_rows, n_cols = set_diameter(n, all_squares)

            SIZE = set_figsize(n_rows, n_cols)

            fig, axes = plt.subplots(n_rows, n_cols, figsize=SIZE)
            fig.suptitle(env_path.parts[-1])

            fig.tight_layout()
            fig.subplots_adjust(top=0.88, bottom=0.5)

            axes = np.array(axes)
            ax = axes.flatten()

            n_axes = ax.size
            if n < n_axes:
                for i in range(1, (n_axes - n) + 1):
                    ax[-i].set_visible(False)

            for idx, mol in enumerate(batch):
                plot_molecule(mol, stats[mol], midplane_cut, water_cut, X, ax[idx])

            mols = [mol.split('_')[1] for mol in list(stats.keys())]
            numbers = '-'.join(mols)
            name = Path(f'{env_path}/{env_path.parts[-1]}_mols-{numbers}.png')
            plt.savefig(name)
    else:
        n = n_mols
        n_rows, n_cols = set_diameter(n, all_squares)

        SIZE = set_figsize(n_rows, n_cols)

        fig, axes = plt.subplots(n_rows, n_cols, figsize=SIZE)
        fig.suptitle(env_path.parts[-1])

        fig.tight_layout()
        fig.subplots_adjust(top=0.9, bottom=0.1)

        axes = np.array(axes)
        ax = axes.flatten()

        n_axes = ax.size
        if n < n_axes:
            for i in range(1, (n_axes - n) + 1):
                ax[-i].set_visible(False)

        for idx, mol in enumerate(list(stats.keys())):
            plot_molecule(mol, stats[mol], midplane_cut, water_cut, X, ax[idx])

        # plt.show()
        # mols = [mol.split('_')[1] for mol in list(stats.keys())]
        # numbers = '-'.join(mols)
        name = Path(f'{env_path}/{env_path.parts[-1]}-{env_path.parts[-2]}.png')
        plt.savefig(name)


def generate_histogram(df):

    stats = dict()

    for idx, elem in enumerate(df.columns.tolist()):
        env = elem[0]
        mol = elem[1]

        X = np.linspace(start=1, stop=10, num=10)

        if env is 'POPG':
            lower = 0.8
            upper = 2.4
        else:
            lower = 1.2
            upper = 2.8

        edges = np.linspace(start=lower, stop=upper, num=11)

        counts, bin_edges = np.histogram(df[env][mol].dropna(), bins=edges, density=True)
        inner = np.fromiter([np.sum(counts[:i]) for i, val in enumerate(counts, 0)], float)
        outer = np.fromiter([np.sum(counts[-i:]) for i, val in enumerate(counts, 1)], float)[::-1]

        f = interp1d(inner, X, assume_sorted=False)
        f_y50 = (inner[-1] - inner[0]) * 0.5
        f_xnew = float(f(f_y50))

        results = [(f_xnew, f_y50), inner, outer]
        stats[mol] = results

    return stats, env, X


def check_positions(stats, env, midplane_cut, water_cut, logfile, environment):

    is_interfacial = False
    for mol, val in stats.items():
        if midplane_cut <= val[0][0] <= water_cut:
            logfile.at[mol, 'interfacial'] = True
            logfile.at[mol, 'y=50%'] = val[0][0]
            is_interfacial = True
        else:
            logfile.at[mol, 'interfacial'] = False
            logfile.at[mol, 'y=50%'] = val[0][0]

    logfile.columns = pd.MultiIndex.from_product([[env], logfile.columns])
    logfile.to_pickle(f'{environment}/{mol}-{env}-log.pickle')

    return is_interfacial


def calculate_probabilities(environment):

    midplane_cut, water_cut = 3.7, 7.94

    environment = environment.resolve()
    mol = environment.parts[-2]
    logfile = pd.DataFrame(columns=['y=50%', 'interfacial'], index=[mol])
    # print(logfile)

    res_dirs = [mol_lambda for mol_lambda in sorted(environment.iterdir()) if mol_lambda.is_dir()]
    df = read_data(res_dirs, mol)
    df.columns = pd.MultiIndex.from_product([[environment.stem], df.columns])
    # print(df)

    stats, env, X = generate_histogram(df)
    plot_molecules(stats, environment, midplane_cut, water_cut, X)
    result = check_positions(stats, env, midplane_cut, water_cut, logfile, environment)

    if result:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-env', '--environment', type=Path, help='Need: path to environment directory with PG '
                                                                 'simulation results')
    # parser.add_argument('-log', '--logfile', type=Path, help='Need: path to pandas dataframe for decision logs.')
    args = parser.parse_args()
    calculate_probabilities(args.environment)
