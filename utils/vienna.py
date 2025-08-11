import os
import sys
import numpy as np
from utils.structure import extract_pairs

# sys.path.append(os.environ.get('VIENNAPATH')) # deprecated due to the python package of ViennaRNA

import RNA


def base_pair_probs(seq, sym=False, scale=True, scale_energy=None):
    fc = RNA.fold_compound(seq)
    if scale:
        if scale_energy is None:
            _, scale_energy = fc.mfe()
        fc.exp_params_rescale(scale_energy)
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    if sym:
        bpp += bpp.T
        unpair = 1 - np.sum(bpp, axis=1)
        bpp[range(len(bpp)), range(len(bpp))] = unpair
    return bpp


def ensemble_defect(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        energy = fc.eval_structure(ss)
        fc.exp_params_rescale(energy)
    fc.pf()
    fc.bpp()
    ed = fc.ensemble_defect(ss)
    return ed


def position_defect(seq, ss, scale=True):
    scale_energy = None
    if scale:
        scale_energy = RNA.fold_compound(seq).eval_structure(ss)
    bpp = base_pair_probs(seq, sym=True, scale=True, scale_energy=scale_energy)
    pairs = extract_pairs(ss)
    defect_pos = [1 - bpp[i, j] for i, j in enumerate(pairs)]
    return defect_pos


def position_defect_mfe(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    ss_mfe_list = subopt(seq)["ss_list"]  # a list of (mfe, structure) tuples
    mfe = ss_mfe_list[0][0]  # minimum free energy
    ss_list = [ss_mfe[1] for ss_mfe in ss_mfe_list]
    if scale:
        fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    sym = True
    if sym:
        bpp += bpp.T
        unpair = 1 - np.sum(bpp, axis=1)
        bpp[range(len(bpp)), range(len(bpp))] = unpair
    pairs = extract_pairs(ss)
    defect_pos = [1 - bpp[i, j] for i, j in enumerate(pairs)]

    return defect_pos, ss_list


def position_ed_pd(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        energy = fc.eval_structure(ss)
        fc.exp_params_rescale(energy)
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    sym = True
    if sym:
        bpp += bpp.T
        unpair = 1 - np.sum(bpp, axis=1)
        bpp[range(len(bpp)), range(len(bpp))] = unpair
    pairs = extract_pairs(ss)
    defect_pos = [1 - bpp[i, j] for i, j in enumerate(pairs)]
    pr = fc.pr_structure(ss)
    pd = -pr
    return defect_pos, pd


def position_ed_pd_mfe(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    ss_mfe_list = subopt(seq)["ss_list"]  # a list of (mfe, structure) tuples
    mfe = ss_mfe_list[0][0]  # minimum free energy
    ss_list = [ss_mfe[1] for ss_mfe in ss_mfe_list]
    if scale:
        fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    sym = True
    if sym:
        bpp += bpp.T
        unpair = 1 - np.sum(bpp, axis=1)
        bpp[range(len(bpp)), range(len(bpp))] = unpair
    pairs = extract_pairs(ss)
    defect_pos = [1 - bpp[i, j] for i, j in enumerate(pairs)]
    pr = fc.pr_structure(ss)
    pd = -pr
    return defect_pos, pd, ss_list


def position_ed_ned_mfe(seq, ss):
    defect_list, ss_list = position_defect_mfe(seq, ss)
    ned = sum(defect_list) / len(defect_list)
    return defect_list, ned, ss_list


def energy(seq, ss):
    fc = RNA.fold_compound(seq)
    return fc.eval_structure(ss)


def mfe(seq):
    fc = RNA.fold_compound(seq)
    ss = fc.mfe()
    return ss


def prob(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(ss)
    return pr


def prob_defect(seq, ss):
    return 1 - prob(seq, ss)


# Print a subopt result as FASTA record
def print_subopt_result(structure, energy, data):
    ss_list = []
    if not structure == None:
        # print(">subopt {:d}".format(data['counter']))
        # print("{}\n{} [{:6.2f}]".format(data['sequence'], structure, energy))
        data["ss_list"].append((energy, structure))
        # increase structure counter
        data["counter"] = data["counter"] + 1


def subopt(seq, e=0):
    subopt_data = {"counter": 0, "sequence": seq, "ss_list": []}
    fc = RNA.fold_compound(seq)
    fc.subopt_cb(e, print_subopt_result, subopt_data)
    subopt_data["ss_list"] = sorted(subopt_data["ss_list"])
    return subopt_data


if __name__ == "__main__":
    if len(sys.argv) > 1:
        rna = sys.argv[1]
    else:
        rna = "AAAAAAAAAACCGCAAAAGCGGGGCCUAAUGGCCGCGGAAUCCGC"
    ss_mfe = mfe(rna)[0]
    bpp = base_pair_probs(rna, sym=True)
    defect = ensemble_defect(rna, ss_mfe)
    defect_pos = position_defect(rna, ss_mfe)
    pr = prob(rna, ss_mfe)
    subopt_data = subopt(rna)
    defect_pos_2, pr_2 = position_ed_pd(rna, ss_mfe)
    _, ned, _ = position_ed_ned_mfe(rna, ss_mfe)
    print("rna:", rna)
    print("mfe:", ss_mfe)
    print("prb:", pr)
    print("pr2:", pr_2)
    print("bpp:", bpp)
    print("Vinenna NED:", defect)
    print("Scratch NED:", sum(defect_pos) / len(defect_pos))
    print("Position ED:", sum(defect_pos_2) / len(defect_pos_2))
    print("        NED:", ned)
    print("subopt:", subopt_data)
