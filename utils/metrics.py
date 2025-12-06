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


def read_structures_from_txt(path):
    structure_list = []
    with open(path, "r") as f:
        for line in f:
            structure_list.append(line.strip())
    return structure_list


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
    filtered_structures = None
    if args.path_filter:
        filtered_structures = read_structures_from_txt(args.path_filter)
    data = []
    matrix_mfe = np.zeros((num_puzzles, len(df_list)), dtype=int)
    matrix_umfe = np.zeros((num_puzzles, len(df_list)), dtype=int)
    for i in range(num_puzzles):
        structure = df_list[0].structure.iloc[i]
        if "(..)" in structure or "(.)" in structure or "()" in structure:  # skip structures with sharp turns
            continue
        if filtered_structures and structure in filtered_structures:
            continue
        count_soved_by_mfe = 0
        count_soved_by_umfe = 0
        objective_list = []
        dist_best_list = []
        ned_best_list = []
        time_list = []
        for j, df_one in enumerate(df_list):
            mfe_list = eval(df_one.mfe_list.iloc[i])
            umfe_list = eval(df_one.umfe_list.iloc[i])
            if mfe_list:
                count_soved_by_mfe += 1
                matrix_mfe[i, j] = 1
            if umfe_list:
                count_soved_by_umfe += 1
                matrix_umfe[i, j] = 1
            objective_list.append((df_one.objective.iloc[i]))
            dist_best_list.append(eval(df_one.dist_best.iloc[i]))
            ned_best_list.append(eval(df_one.ned_best.iloc[i]))
            time_list.append(df_one.time.iloc[i])
        dist_best, _ = min(dist_best_list)
        ned_best, _ = min(ned_best_list)
        time_mean = np.mean(time_list)
        if dist_best < 0:
            dist_best = 0
        data.append(
            [
                i,
                structure,
                min(objective_list),
                np.mean(objective_list),
                count_soved_by_mfe,
                count_soved_by_umfe,
                count_soved_by_mfe > 0,
                count_soved_by_umfe > 0,
                dist_best,
                ned_best,
                time_mean,
            ]
        )
    df_joint = pd.DataFrame(
        data,
        columns=(
            "index",
            "structure",
            "objective",
            "objective_mean",
            "count_solved_mfe",
            "count_solved_umfe",
            "is_solved_mfe",
            "is_solved_umfe",
            "dist_best",
            "ned_best",
            "time_mean",
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
    print()

    print("best distance:")
    print("-------------------------------------")
    print(f"mean: {df_joint.dist_best.mean():.2f} std: {df_joint.dist_best.std():.2f}")
    print()

    print("best ned:")
    print("-------------------------------------")
    print(f"mean: {df_joint.ned_best.mean():.4f} std: {df_joint.ned_best.std():.4f}")
    print()

    print("time statistics:")
    print("-------------------------------------")
    print(f"time mean: {df_joint.time_mean.mean():.2f} std: {df_joint.time_mean.std():.2f}")
    print()

    print("Objective statistics:")
    print("-------------------------------------")
    print(f"objective arithmic mean: {df_joint.objective.mean():.4f}")
    print(f"1 - obj.  arithmic mean: {1 - df_joint.objective.mean():.4f}")
    # print(f"objective arithmic std: {df_joint.objective.std():.4f}")
    print()

    objective_list = df_joint.objective
    # geometric mean and std
    geometric_mean = np.exp(np.log(objective_list).mean())
    # geometric_std = np.exp(np.log(prob_list).std())
    if geometric_mean > 1e-4:
        print(f"objective geometric mean: {geometric_mean:.4f}")
    else:  # scientific notation for very small mean value
        print(f"objective geometric mean: {geometric_mean:.4e}")

    objective_complement_list = 1 - df_joint.objective
    geometric_mean_complement = np.exp(np.log(objective_complement_list).mean())
    if geometric_mean_complement > 1e-4:
        print(f"1 - obj.  geometric mean: {geometric_mean_complement:.4f}")
    else:  # scientific notation for very small mean value
        print(f"1 - obj.  geometric mean: {geometric_mean_complement:.4e}")
    print()

    print("time statistics:")
    print("-------------------------------------")
    print(f"average: {df_joint.time_mean.mean():.2f} std: {df_joint.time_mean.std():.2f}")
    print()

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
        obj_best_list = []
        for j, df_one in enumerate(df_list):
            dist_best = eval(df_one.dist_best.iloc[i])
            dist_best_list.append(dist_best)
            obj_best = (df_one.objective.iloc[i], df_one.rna.iloc[i])  # actually the first element is not a distance but an objective value
            obj_best_list.append(obj_best)
        obj_best = min(obj_best_list)
        dist_best = min(dist_best_list)
        dist, _ = dist_best 
        _, x = obj_best
        if dist >= 0:
            y_mfe, d_best = argmin_dist(x, structure)
            # assert d_best == dist, f"Distance mismatch: {d_best} != {dist}"
            assert d_best >= dist, f"Distance mismatch: {d_best} < {dist}"
            data.append([i, structure, d_best, dist, x, y_mfe])
        else:
            data.append([i, structure, dist, dist, x, structure])  # no valid sequence found
    df_joint = pd.DataFrame(data, columns=("index", "y", "dist", "dist_raw", "x", "y_mfe"))
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
    parser.add_argument("--path_filter", type=str, default="", help="path to a txt file containing structures to be filtered out")

    args = parser.parse_args()
    print(f"args: {args}")

    main(args)
