"""This script reads multiple CSV files containing MFE (Minimum Free Energy) data,
counts the number of MFE and UMFE (Unconstrained Minimum Free Energy) solutions,
and prints various statistics about the solutions.
Example usage:
python metrics.py --path "samfeo_output*.csv" or python metrics.py --path samfeo_output\*.csv
"""

import argparse
import glob

import numpy as np
import pandas as pd


def read_dataframes(path):
    df_list = []
    print()
    print("Reading files ..")
    print("-------------------------------------")
    for filename in glob.glob(path):
        print(filename, end=" ")
        df_onerun = pd.read_csv(filename)
        df_list.append(df_onerun)
        print(len(df_onerun))
    print("-------------------------------------")
    print("Done reading files.")
    return df_list


def count_mfe(df_list):
    assert all([len(df) == len(df_list[0]) for df in df_list])
    num_puzzles = len(df_list[0])
    print(f"num_puzzles: {num_puzzles}")
    data = []
    matrix_mfe = np.zeros((num_puzzles, len(df_list)), dtype=int)
    matrix_umfe = np.zeros((num_puzzles, len(df_list)), dtype=int)
    for i in range(num_puzzles):
        structure = df_list[0].structure.iloc[i]
        count_mfe = 0
        count_umfe = 0
        for j, df_one in enumerate(df_list):
            mfe_list = eval(df_one.mfe_list.iloc[i])
            umfe_list = eval(df_one.umfe_list.iloc[i])
            if mfe_list:
                count_mfe += 1
                matrix_mfe[i, j] = 1
            if umfe_list:
                count_umfe += 1
                matrix_umfe[i, j] = 1

        data.append(
            [i, structure, count_mfe, count_umfe, count_mfe > 0, count_umfe > 0]
        )
    df_joint = pd.DataFrame(
        data,
        columns=(
            "index",
            "structure",
            "count_mfe",
            "count_umfe",
            "has_mfe",
            "has_umfe",
        ),
    )

    print()
    print("MFE metrics:")
    print("-------------------------------------")
    print(f"solved by  mfe: {df_joint.has_mfe.sum()}")
    print(f"solved by umfe: {df_joint.has_umfe.sum()}")
    print(
        f"  mfe mean and std: {matrix_mfe.sum(axis=0).mean():.1f} {matrix_mfe.sum(axis=0).std():.1f}"
    )
    print(
        f" umfe mean and std: {matrix_umfe.sum(axis=0).mean():.1f} {matrix_umfe.sum(axis=0).std():.1f}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default="")

    args = parser.parse_args()
    print(f"args: {args}")

    df_list = read_dataframes(args.path)
    count_mfe(df_list)
