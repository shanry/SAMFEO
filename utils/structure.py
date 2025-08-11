"""
This module provides functions to handle RNA secondary structures,
including extracting pairs, calculating distances, etc.
"""


def extract_pairs(ss):
    """Extract pairs from a secondary structure string.
    Returns a list of indices where each position has the index of its pair, or itself if unpaired.
    """
    pairs = list(range(len(ss)))
    stack = []
    for i, c in enumerate(ss):
        if c == ".":
            pass
        elif c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs


def extract_pairs_list(ss):
    """Extract pairs from a secondary structure string.
    Returns a list of tuples where each tuple contains the indices of paired positions.
    """
    pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c == ".":
            pass
        elif c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.append((j, i))
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs


def pairs_match(ss):
    """Extract pairs in a secondary structure string.
    Returns a dictionary where keys are indices of left/right brackets and values are indices of right/left brackets.
    """
    assert len(ss) > 5
    pairs = dict()
    stack = []
    for i, s in enumerate(ss):
        if s == ".":
            pass
        elif s == "(":
            stack.append(i)
        elif s == ")":
            j = stack.pop()
            assert j < i
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(
                f"the value of structure at position: {i} is not right: {s}!"
            )
    return pairs


def struct_dist(s1, s2):
    assert len(s1) == len(s2), f"len(s1)={len(s1)}, len(s2)={len(s2)}"
    pairs_1 = pairs_match(s1)
    pairs_2 = pairs_match(s2)
    union = len(pairs_1.keys() | pairs_2.keys())
    overlap = len(s1) - union
    for k in pairs_1:
        if k in pairs_2 and pairs_1[k] == pairs_2[k]:
            overlap += 1
    return len(s1) - overlap


if __name__ == "__main__":
    ss = "..........((((....))))((((....))))((((...))))"
    pairs = extract_pairs(ss)
    print("structure:", ss)
    print("pairs:", pairs)
    pairs_list = extract_pairs_list(ss)
    print("pairs_list: ", pairs_list)
    pairs_dict = pairs_match(ss)
    print("pairs_dict: ", pairs_dict)
