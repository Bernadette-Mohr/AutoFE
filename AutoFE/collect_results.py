import pandas as pd
from pathlib import Path
import argparse


def create_df(directory, dataframe, _round):
    _dir = directory.resolve()
    first_step = _dir / 'POPG'
    delta = u'\u0394'

    if not _dir.is_dir():
        print(f'{_dir} does not exist!')
        return
    elif not first_step.is_dir():
        print(f'{first_step} does not exist!')
        return
    else:
        if dataframe is None and _round is not None:
            molecules = [mol.stem for mol in sorted(first_step.iterdir()) if mol.is_dir()]
            df = pd.DataFrame(index=molecules, columns=['y=50%', 'interfacial', 'OCTANE', 'POPG', 'WATER', 'PG O->I',
                                                        'PG W->I', 'CDL2', 'CL O->I', 'CL W->I', f'{delta}{delta}G PG->CL'])
            results = _dir / f'results-{_round}.pickle'
        else:
            df = pd.read_pickle(dataframe)
            results = dataframe

        for mol in first_step.iterdir():
            if mol.is_dir():
                idx = mol.stem
                first_df = pd.read_pickle(mol / f'{idx}-POPG-log.pickle')
                hist_results = first_df.columns.get_level_values(1).to_list()
                if hist_results[0] not in df and hist_results[1] not in df:
                    df.at[idx, hist_results[0]] = first_df.at[idx, ('POPG', hist_results[0])]
                    df.at[idx, hist_results[1]] = first_df.at[idx, ('POPG', hist_results[1])]
                else:
                    if pd.isna(df.loc[idx, hist_results[0]]):
                        df.at[idx, hist_results[0]] = first_df.at[idx, ('POPG', hist_results[0])]
                        df.at[idx, hist_results[1]] = first_df.at[idx, ('POPG', hist_results[1])]

        archives = [pkl for pkl in sorted(_dir.glob('*.pickle')) if pkl.is_file()]
        for archive in archives:
            content = pd.read_pickle(archive)
            idx = content.index.to_list()[0]
            for column in content.columns:
                if column not in df:
                    df.at[idx, column] = content.at[idx, column]
                else:
                    if pd.isna(df.loc[idx, column]):
                        df.at[idx, column] = content.at[idx, column]

        df.to_pickle(results)
        df.to_csv(_dir / f'results-{_round}.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect all individual FEs of a single learning round into one '
                                                 'pandas dataframe.')
    parser.add_argument('-dir', '--directory', type=Path, required=True,
                        help='Directory with current round of simulations.')
    parser.add_argument('-df', '--dataframe', type=Path, required=False, default=None,
                        help='If results finished later and have to be added, pass the existing pandas dataframe.')
    parser.add_argument('-r', '--round', type=str, required=False, default=None,
                        help='name of dataframe to be created, eg. \"round_0\".')

    args = parser.parse_args()

    create_df(args.directory, args.dataframe, args.round)
