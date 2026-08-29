# SAMFEO

Implementation of the RNA design method from the paper (ISMB 2023):

[1] Zhou, T., Dai, N., Li, S., Ward, M., Mathews, D.H. and Huang, L., 2023. RNA design via structure-aware multifrontier ensemble optimization. Bioinformatics, 39(Supplement_1), pp.i563-i571.

The main branch is under active development and is continuously updated to improve numerical precision, code readability, and functionality. The reference implementation used for ISMB 2023 is preserved in the `ismb2023` branch: https://github.com/shanry/SAMFEO/tree/ismb2023

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/shanry/SAMFEO
    cd SAMFEO
    ```
2.  **Install dependencies:**
    This project is developed using Python 3.11, and should also work with other Python versions.
    ```bash
    pip install -r requirements.txt
    ```

## Testing

```bash
pytest main.py -s
```

## Quick Start

### Online mode (input a single structure at a time)

```bash
python main.py --online  # plus other design options
```

Example:
```
$ echo "(((((......)))))" | python main.py --online --step 100
args:
Namespace(path='', object='pd', k=10, t=1, step=100, stay=2000, name='', init='cg', repeat=1, start=0, nomfe=False, nosm=False, log=False, online=True, worker_count=1, batch_size=1)
(((((......)))))
seed_np: 2020
steps: 100, t: 1, k: 10, structured mutation: True, objective name: probability_defect
name_pair: cg
pair_pool: ['CG', 'GC']
0 CGCGCAAAAAAGCGCG:  8.0927e-01
1 CCGCGAAAAAACGCGG:  9.1725e-01
2 CCCGGAAAAAACCGGG:  9.4103e-01
3 GGCCGAAAAAACGGCC:  9.4508e-01
4 CCCCCAAAAAAGGGGG:  9.3840e-01
5 GCGGCAAAAAAGCCGC:  9.7389e-01
6 GCCCCAAAAAAGGGGC:  9.7457e-01
7 GCGCCAAAAAAGGCGC:  9.7389e-01
8 GGGGCAAAAAAGCCCC:  9.7455e-01
9 GCCCGAAAAAACGGGC:  9.4502e-01
iter:    10      value: -9.7457e-01      mfe count:    18        umfe count: 18  best iter: 0 improve: 0.00e+00
iter:    20      value: -9.7457e-01      mfe count:    26        umfe count: 26  best iter: 0 improve: 0.00e+00
iter:    30      value: -9.7585e-01      mfe count:    35        umfe count: 35  best iter: 28 improve: -1.28e-03
iter:    40      value: -9.7585e-01      mfe count:    42        umfe count: 42  best iter: 28 improve: 0.00e+00
iter:    50      value: -9.7585e-01      mfe count:    51        umfe count: 51  best iter: 28 improve: 0.00e+00
iter:    60      value: -9.7585e-01      mfe count:    59        umfe count: 59  best iter: 28 improve: 0.00e+00
iter:    70      value: -9.7585e-01      mfe count:    68        umfe count: 68  best iter: 28 improve: 0.00e+00
iter:    80      value: -9.7585e-01      mfe count:    77        umfe count: 77  best iter: 28 improve: 0.00e+00
iter:    90      value: -9.7585e-01      mfe count:    85        umfe count: 85  best iter: 28 improve: 0.00e+00
iter:   100      value: -9.7585e-01      mfe count:    92        umfe count: 92  best iter: 28 improve: 0.00e+00
RNA sequence: 
GCCCCGAAAAAGGGGC
ensemble objective:  0.02414993696568135
(((((......)))))
(((((......)))))
structure distance: 0
count of mfe solutions: 92
count of umfe solutions: 92
[RNAStructure('GGGGCAAAAACGCCCC', 0.9745324428444898), RNAStructure('GGGGCAGAACAGCCCC', 0.9745537982728839), RNAStructure('GGGGCAAAGAAGCCCC', 0.9745522394520901), RNAStructure('GGGGCAGAAAAGCCCC', 0.9745541635522057), RNAStructure('GGCCCAAAAAAGGGCC', 0.974676524461917), RNAStructure('GCCCCAAUAAAGGGGC', 0.9745688330920909), RNAStructure('GCCCCAAAAAAGGGGC', 0.9745692982930593), RNAStructure('GCCCCGAAAAAGGGGC', 0.9758500630343186), RNAStructure('GGGGCAAAAAAGCCCC', 0.9745543077208556), RNAStructure('GGAGCGAAGAAGCUCC', 0.9756107700034071)]
 mfe samples: ['GGGUGCAAGACCACCC', 'GGGGCAAAGAAGCUUC', 'GGGGCAGAACAGCCCC', 'GGUGCAGAAAAGUGCC', 'GCCCCAAAAGAGGGGC']

umfe samples: ['GGGUGCAAGACCACCC', 'GGGGCAAAGAAGCUUC', 'GGGGCAGAACAGCCCC', 'GGUGCAGAAAAGUGCC', 'GCCCCAAAAGAGGGGC']

best: GCCCCGAAAAAGGGGC:  9.7585e-01

ned_best: (0.0030683500152726695, 'GCCCCGAAAAAGGGGC')

prob_best: (0.9758500630343186, 'GCCCCGAAAAAGGGGC')

dist_best: (0, 'GCCCCAAAAGAGGGGC')

full results are saved in the file: results_20250718161920_1975123.json  # results_{timestamp}_{random_id}.json
```

### Batch mode (input a file with one structure per line)

```bash
python main.py --t 1 --k 10 --object pd --path data/eterna/eterna100.txt  # add --batch_size 50 --worker_count 20 for multiprocessing speedups
```

### Replicate reported results (5 repeats)

```bash
python main.py --t 1 --k 10 --object pd --path data/eterna/eterna100.txt --repeat 5  # add --batch_size 50 --worker_count 20 for multiprocessing speedups
```

### Design options

| Option | Description |
| --- | --- |
| `--online` | Enable online mode (read a single structure from stdin). |
| `--t` | Sampling temperature T. Higher T increases randomness; lower T is greedier. Default: `1`. |
| `--k` | Frontier (priority queue) size. Larger values explore more candidates per iteration at higher cost. Default: `10`. |
| `--object` | Optimization objective: `pd` (probability defect) or `ned` (normalized ensemble defect). Default: `pd`. |
| `--init` | Sequence initialization method: `cg` (constraint-guided/targeted) or `all` (uniform random). Default: `cg`. |
| `--nosm` | Disable structured mutation (ablates structure-aware grouped changes). |
| `--nomfe` | Disable "MFE as product"; only the single best-scoring sequence is kept as the designed output. |
| `--step` | Maximum number of iterations. Default: `5000`. |
| `--repeat` | Number of repeated runs (batch mode only). Default: `1`. |
| `--start` | Starting index for naming output files (batch mode only). Default: `0`. |
| `--batch_size` | Number of structures processed concurrently in batch mode. |
| `--worker_count` | Number of CPU workers for batch mode; must be ≤ `--batch_size`. |

### Batch mode output format

The design results are written to a CSV file with the following columns:

| Column | Description |
| --- | --- |
| `structure` | The input puzzle (target secondary structure). |
| `objective` | The best (lowest) objective value achieved during optimization. |
| `rna` | The sequence with the best objective during optimization. |
| `mfe` | The MFE structure of `rna`, given by the ViennaRNA folding engine. |
| `dist` | The distance between `mfe` and `structure`. |
| `mfe_list` | A list of all the MFE solutions found, ordered by the time each was detected. |
| `umfe_list` | A list of all the uMFE solutions found, ordered by the time each was detected. |
| `k_best` | The priority queue at the final iteration. |
| `ned_best` | The best NED value and the corresponding sequence. |
| `prob_best` | The best probability value and the corresponding sequence. |
| `dist_best` | The best distance value and the corresponding sequence. |
| `log` | The objective values at each iteration. |
| `time` | The total time used to design the input `structure`. |
