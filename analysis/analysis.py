import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.preprocessing.subsampling import statistical_inefficiency
from pymbar.mbar import MBAR


def skip_equilibration_time(df, skiptime):
    _skip = df.index.get_level_values(0).get_loc(skiptime)
    # print(_skip)
    return df.iloc[_skip:]


def read_potentials(res_dirs, skiptime, TEMP):
    rgx = re.compile(r'prod-\d{2,3}\.xvg')
    xvg_files = list()

    for res_dir in res_dirs:
        [xvg_files.append(_file) for _file in res_dir.iterdir() if re.match(rgx, _file.name)]

    if not xvg_files:
        return
    else:
        print(f'\nGathering u_nk from {len(xvg_files)} files')
        u_nk = [extract_u_nk(_file, T=TEMP) for _file in xvg_files]

        print(f'\nRemoving first {skiptime} ns as equilibration time')
        u_nk = [skip_equilibration_time(df, skiptime) for df in u_nk]
        print(f'\nSubsampling u_nk on statistical inefficiency')
        u_nk = pd.concat([statistical_inefficiency(df, column=df.iloc[:, 0]) for df in u_nk])

        return u_nk


"""
MBAR returns energies in kBT, have to be converted to kcal/mol.
"""


def convert_FEunits(Deltaf_ij, dDeltaf_ij, K, TEMP):
    kB = 1.3806488 * 6.02214129 / 1000.0  # Boltzmann's constant (kJ/mol/K).
    beta = 1. / (kB * TEMP)
    beta_report = 4.184 * beta

    return Deltaf_ij[0:K - 1] / beta_report, dDeltaf_ij[0:K - 1] / beta_report


def perform_analysis(data, TEMP):
    df_allk, ddf_allk = list(), list()
    N_k = list()

    if len(data.index.names) == 2:
        # time = data.index.get_level_values(0)
        vdw_lambda = data.index.get_level_values(1)
        n_vdw = len(vdw_lambda.unique())
        lambdas = sorted(set(vdw_lambda))
        for _lambda in lambdas:
            N_k.append(len(data.xs(_lambda, level='vdw-lambda')))

    else:
        # time = data.index.get_level_values(0)
        coul_lambda = data.index.get_level_values(1)
        vdw_lambda = data.index.get_level_values(2)
        # n_coul = len(coul_lambda.unique())
        n_vdw = len(vdw_lambda.unique())
        lambda_pairs = sorted(set(tuple([m, n]) for m, n in zip(coul_lambda, vdw_lambda)))
        for pair in lambda_pairs:
            N_k.append(len(data.xs(pair, level=('coul-lambda', 'vdw-lambda'))))

    u_kn = data.to_numpy().T
    # print(u_kn.shape)

    # print(len(N_k), sum(N_k))
    K = len(N_k)

    print(f'\nEstimating the free energy change with MBAR...')
    mbar = MBAR(u_kn, N_k, verbose=True, relative_tolerance=1e-10, initialize='BAR')
    Deltaf_ij, dDeltaf_ij, theta_ij = mbar.getFreeEnergyDifferences(uncertainty_method='svd-ew', return_theta=True)

    print(f'\nConvert dimensionless free energies to kcal/mol...')
    Deltaf_ij, dDeltaf_ij = convert_FEunits(Deltaf_ij, dDeltaf_ij, K, TEMP)

    print(f'\nEstimating DeltaG and errors for all pairs of adjacent states...')
    for k in range(len(N_k) - 1):
        df_allk = np.append(df_allk, Deltaf_ij[k, k + 1])
        ddf_allk = np.append(ddf_allk, np.nan_to_num(dDeltaf_ij[k, k + 1]))

    return df_allk, ddf_allk, n_vdw, K


def print_logfile(env, dF, ddF, skiptime, TEMP):
    dF_list = list(np.round(dF['TOTAL'], 3))
    ddF_list = list(np.round(ddF['TOTAL'], 3))
    results = zip(dF_list, ddF_list)
    logpath = env / f'MBAR-results_{env.resolve().parts[-1]}-{env.resolve().parts[-2]}.txt'

    with open(logpath, 'a+') as logfile:
        logfile.write(f'Free energy estimates done using alchemlyb and pymbar. Simulations run at T = {TEMP} K. '
                      f'First {skiptime} ns were removed as equilibration time.\n\n')
        logfile.write(f'    step    \t  MBAR [kcal/mol]  \n')
        logfile.write(f'---------------------------------\n')
        for i, (val, err) in enumerate(results):
            logfile.write(f'{i + 1:3} -> {i + 2:3}\t{val:6} +- {err:5}\n')
        logfile.write(f'---------------------------------\n')
        logfile.write(f'Coulomb:\t{np.round(np.sum(dF["Coulomb"]), 3):>7} +- '
                      f'{np.round(np.power(np.sum(ddF["Coulomb"]), 2), 3)}\n')
        logfile.write(f'VdWaals:\t{np.round(np.sum(dF["VdWaals"]), 3):>7} +- '
                      f'{np.round(np.power(np.sum(ddF["VdWaals"]), 2), 3)}\n')
        logfile.write(f'  TOTAL:\t{np.round(np.sum(dF["TOTAL"]), 3):>7} +- '
                      f'{np.round(np.power(np.sum(ddF["TOTAL"]), 2) * 0.5, 3)}\n')


def calculate_FEs(env, df, ddf, n_vdw, K, skiptime, TEMP, env_mol=None):
    startcoul = n_vdw - 1
    endcoul = len(df)
    startvdw = 0
    endvdw = n_vdw - 1

    segments = ['Coulomb', 'VdWaals', 'TOTAL']
    segmentstarts = [startcoul, startvdw, 0]
    segmentends = [endcoul, endvdw, K - 1]

    dF, ddF = dict(), dict()

    for seg in range(len(segments)):
        segment = segments[seg]
        segstart = segmentstarts[seg]
        segend = segmentends[seg]

        dF[segment] = df[segstart:segend]
        ddF[segment] = ddf[segstart:segend]

    total = np.round(np.sum(dF["TOTAL"]), 3)
    d_total = np.round(np.power(np.sum(ddF["TOTAL"]), 2) * 0.5, 3)

    if env_mol is not None:
        print_logfile(env_mol, dF, ddF, skiptime, TEMP)
    else:
        print_logfile(env, dF, ddF, skiptime, TEMP)

    return total, d_total


def analyze_FE(environments, skiptime, archive, molecules):
    """
    Martini is parameterized at T = 300 K, so hard coding the simulation temperature is justifiable since it will
    always be 300 K.
    """
    TEMP = 300

    path = environments[0].resolve().parent

    if not molecules:

        mol = environments[0].resolve().parts[-2]

        if not archive:
            dataframe = pd.DataFrame(columns=[env.stem for env in sorted(environments) if env.is_dir()], index=[mol])
        else:
            dataframe = pd.read_pickle(archive)

        for env in sorted(environments):
            if env not in dataframe:
                dataframe[env.stem] = ''
            res_dirs = [subdir for subdir in sorted(env.iterdir()) if subdir.is_dir()]
            u_nk = read_potentials(res_dirs, skiptime, TEMP)
            if u_nk is not None:
                df, ddf, n_vdw, K = perform_analysis(u_nk, TEMP)
                total, d_total = calculate_FEs(env, df, ddf, n_vdw, K, skiptime, TEMP)
                dataframe.at[mol, env.stem] = (total, d_total)

        dataframe.to_pickle(f'{path}/FEs-{mol}.pickle')

    else:
        for molecule in sorted(molecules.iterdir()):
            mol = molecule.resolve().stem

            archive = path / f'FEs-{mol}.pickle'
            if archive.is_file():
                dataframe = pd.read_pickle(archive)
            else:
                dataframe = pd.DataFrame(columns=[env.stem for env in sorted(environments) if env.is_dir()],
                                         index=[mol])

            for env in sorted(environments):
                env_mol = env / mol
                if env not in dataframe:
                    dataframe[env.stem] = ''
                res_dirs = [subdir for subdir in sorted(env_mol.iterdir()) if subdir.is_dir()]
                u_nk = read_potentials(res_dirs, skiptime, TEMP)
                if u_nk is not None:
                    df, ddf, n_vdw, K = perform_analysis(u_nk, TEMP)
                    total, d_total = calculate_FEs(env, df, ddf, n_vdw, K, skiptime, TEMP, env_mol)
                    dataframe.at[mol, env.stem] = (total, d_total)

            dataframe.to_pickle(f'{path}/FEs-{mol}.pickle')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Automatic MBAR analysis of free energies.')
    parser.add_argument('-env', '--environment', type=Path, nargs='+', required=True,
                        help='Need: Environment directories in which FE calculations were run.')
    parser.add_argument('-s', '--skiptime', type=float,
                        help='Discard data prior to this specified time as \'equilibration\' data. Units picoseconds. '
                             'Default: 0 ps.', default=0)
    parser.add_argument('-a', '--archive', type=Path, required=False,
                        help='If it has been created in a previous step, pass path to FE results archive (*.pickle)',
                        default=None)
    parser.add_argument('-mol', '--molecules', type=Path, required=False,
                        help='If FE calculation has to be done in $HOME on the cluster, provide location of molecule '
                             'topologies for iterating.', default=None)

    args = parser.parse_args()

    analyze_FE(args.environment, args.skiptime, args.archive, args.molecules)
