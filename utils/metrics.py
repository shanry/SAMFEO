"""This script reads multiple CSV files output from SAMFEO,
counts the number of solved puzzles / structures by MFE and uMFE (unique Minimum Free Energy) criteria,
and prints various statistics about the solutions.
Example usage:
python metrics.py --path "samfeo_output*.csv" --mfe
or python metrics.py --path samfeo_output\*.csv --mfe

find best distance for each puzzle, along with the corresponding sequence and MFE structure.
Example usage:
python metrics.py --path "samfeo_output*.csv" --dist
or python metrics.py --path samfeo_output\*.csv --dist
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
        count_soved_by_mfe = 0
        count_soved_by_umfe = 0
        objective_list = []
        dist_best_list = []
        for j, df_one in enumerate(df_list):
            mfe_list = eval(df_one.mfe_list.iloc[i])
            umfe_list = eval(df_one.umfe_list.iloc[i])
            if mfe_list:
                count_soved_by_mfe += 1
                matrix_mfe[i, j] = 1
            if umfe_list:
                count_soved_by_umfe += 1
                matrix_umfe[i, j] = 1
            objective_list.append(df_one.objective.iloc[i])
            dist_best_list.append(df_one.dist_best.iloc[i])
        data.append(
            [
                i,
                structure,
                min(objective_list),
                count_soved_by_mfe,
                count_soved_by_umfe,
                count_soved_by_mfe > 0,
                count_soved_by_umfe > 0,
            ]
        )
    df_joint = pd.DataFrame(
        data,
        columns=(
            "index",
            "structure",
            "objective",
            "count_solved_mfe",
            "count_solved_umfe",
            "is_solved_mfe",
            "is_solved_umfe",
        ),
    )

    print()
    print("MFE metrics:")
    print("-------------------------------------")
    print(f"solved by  mfe: {df_joint.is_solved_mfe.sum()}")
    print(f"solved by umfe: {df_joint.is_solved_umfe.sum()}")
    print(
        f"  mfe mean and std: {matrix_mfe.sum(axis=0).mean():.1f} {matrix_mfe.sum(axis=0).std():.1f}"
    )
    print(
        f" umfe mean and std: {matrix_umfe.sum(axis=0).mean():.1f} {matrix_umfe.sum(axis=0).std():.1f}"
    )
    print("Objective statistics:")
    print("-------------------------------------")
    print(f"objective mean: {df_joint.objective.mean():.2f}")

    # save the joint dataframe to a CSV file
    df_joint.to_csv("mfe_counts.csv", index=False)
    print("Joint metrics saved to mfe_counts.csv")


def find_best_distance(df_list):
    def argmin_dist(x, y):
        from utils.vienna import subopt
        from utils.structure import struct_dist

        y_mfe_list = subopt(x)["ss_list"]
        y_best = y_mfe_list[0][1]
        d_best = struct_dist(y, y_best)
        for e, y_mfe in y_mfe_list[1:]:
            d = struct_dist(y, y_mfe)
            if d < d_best:
                d_best = d
                y_best = y_mfe
        return y_best, d_best

    assert all([len(df) == len(df_list[0]) for df in df_list])
    num_puzzles = len(df_list[0])
    print(f"num_puzzles: {num_puzzles}")
    data = []
    for i in range(num_puzzles):
        structure = df_list[0].structure.iloc[i]
        dist_best_list = []
        for j, df_one in enumerate(df_list):
            dist_best = eval(df_one.dist_best.iloc[i])
            dist_best_list.append(dist_best)
        dist_best = min(dist_best_list)
        dist, x = dist_best
        if dist >= 0:
            y_mfe, d_best = argmin_dist(x, structure)
            assert d_best == dist, f"Distance mismatch: {d_best} != {dist}"
            data.append([i, structure, dist, x, y_mfe])
    df_joint = pd.DataFrame(data, columns=("index", "y", "dist", "x", "y_mfe"))
    # save the joint dataframe to a CSV file
    df_joint.to_csv("best_distance.csv", index=False)
    print("Best distances saved to best_distance.csv")


def main(args):
    if args.mfe:
        print("Counting MFE statistics ..")
        df_list = read_dataframes(args.path)
        count_mfe(df_list)
    if args.dist:
        print("Finding best distance ..")
        df_list = read_dataframes(args.path)
        find_best_distance(df_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default="")
    parser.add_argument("--dist", action="store_true", help="get best distance")
    parser.add_argument("--mfe", action="store_true", help="get mfe statistics")

    args = parser.parse_args()
    print(f"args: {args}")

    main(args)
