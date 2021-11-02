from pathlib import Path
import pandas as pd
import sys


def get_difference(df, idx):
    delta = u'\u0394'
    # arrow = '->'

    if 'CDL2' in df and not df['CDL2'].isnull().values.any():
        if 'PG W->I' not in df:
            df[f'PG O->I'] = df.at[idx, 'POPG'][0] - df.at[idx, 'OCTANE'][0]
            df[f'PG W->I'] = df.at[idx, 'POPG'][0] - df.at[idx, 'WATER'][0]
            df[f'CL O->I'] = df.at[idx, 'CDL2'][0] - df.at[idx, 'OCTANE'][0]
            df[f'CL W->I'] = df.at[idx, 'CDL2'][0] - df.at[idx, 'WATER'][0]
            df[f'{delta}{delta}G PG->CL'] = df[f'CL W->I'] - df[f'PG W->I']
        else:
            df[f'CL O->I'] = df.at[idx, 'CDL2'][0] - df.at[idx, 'OCTANE'][0]
            df[f'CL W->I'] = df.at[idx, 'CDL2'][0] - df.at[idx, 'WATER'][0]
            df[f'{delta}{delta}G PG->CL'] = df[f'CL W->I'] - df[f'PG W->I']
    else:
        df[f'PG O->I'] = df.at[idx, 'POPG'][0] - df.at[idx, 'OCTANE'][0]
        df[f'PG W->I'] = df.at[idx, 'POPG'][0] - df.at[idx, 'WATER'][0]

    return df


def load_archive(archive):
    path = archive.resolve()
    df = pd.read_pickle(path)
    idx = df.index.to_list()[0]
    cols = df.columns
    for col in cols:
        if isinstance(df.at[idx, col], str):
            df.at[idx, col] = tuple(float(s) for s in df.at[idx, col].strip("()").split(","))

    if 'OCTANE' not in df or 'WATER' not in df:
        return
    elif df['OCTANE'].isnull().values.any() or df['WATER'].isnull().values.any():
        return
    else:
        df = get_difference(df, idx)

        if 'CDL2' in df and not df['CDL2'].isnull().values.any():
            df.to_pickle(path)

        else:
            if (df['PG O->I'] >= 0).any() or (df['PG W->I'] >= 0).any():
                print(f'\nIn PG not interfacial!\n')
                df.to_pickle(path)
                return sys.exit(1)
            else:
                print(f'\nIn PG interfacial!\n')
                df.to_pickle(path)
                return sys.exit(0)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('Reading the different archive files, doing evaluation.')
    parser.add_argument('-a', '--archive', type=Path, required=True,
                        help='Need: Path to pandas dataframe.pickle.')
    args = parser.parse_args()

    load_archive(args.archive)
