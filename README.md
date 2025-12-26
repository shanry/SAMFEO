# SAMFEO
Implementation of the RNA design methd of the paper (ISMB 2023):

[1] Zhou, T., Dai, N., Li, S., Ward, M., Mathews, D.H. and Huang, L., 2023. RNA design via structure-aware multifrontier ensemble optimization. Bioinformatics, 39(Supplement_1), pp.i563-i571.

## Installation
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/shanry/SAMFEO
    cd SAMFEO
    git checkout samfeo_prec
    ```
2.  **Install dependencies:**
    This project is developed using Python 3.11, and should also work for other Python versions.
    ```bash
    pip install -r requirements.txt
    ```

## Testing
```bash
pytest main.py -s
```

## Quick Start
### Online mode (input single structure each time)
``python main.py --online # plus other design options`` 


Example:
```
$ echo "(((((......)))))" | python main.py --online --step 100
args:
Namespace(path='', object='pd', k=10, t=1, step=100, stay=2000, name='', init='cg', repeat=1, start=0, nomfe=False, nosm=False, bp=False, log=False, online=True, para=False, worker_count=10, batch_size=20)
(((((......)))))
seed_np: 2020
steps: 100, t: 1, k: 10, structured mutation: True, ensemble objective: position_ed_pd_mfe
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
count of mfe solutsion: 92
count of umfe solutions: 92
[RNAStructure('GGGGCAAAAACGCCCC', 0.9745324428444898), RNAStructure('GGGGCAGAACAGCCCC', 0.9745537982728839), RNAStructure('GGGGCAAAGAAGCCCC', 0.9745522394520901), RNAStructure('GGGGCAGAAAAGCCCC', 0.9745541635522057), RNAStructure('GGCCCAAAAAAGGGCC', 0.974676524461917), RNAStructure('GCCCCAAUAAAGGGGC', 0.9745688330920909), RNAStructure('GCCCCAAAAAAGGGGC', 0.9745692982930593), RNAStructure('GCCCCGAAAAAGGGGC', 0.9758500630343186), RNAStructure('GGGGCAAAAAAGCCCC', 0.9745543077208556), RNAStructure('GGAGCGAAGAAGCUCC', 0.9756107700034071)]
 mfe samples: ['GGGCUCAAGAUGGCCC', 'GCCCUUAAAACGGGGC', 'GCCCAAAAAACUGGGC', 'GGCACAGAAAAGUGCC', 'GCCCUAAUAACAGGGC', 'GGGUGCAAGACCACCC', 'GGGGCAAAGAAGCUUC', 'GGGGCAGAACAGCCCC', 'GGUGCAGAAAAGUGCC', 'GCCCCAAAAGAGGGGC']

umfe samples: ['GGGCUCAAGAUGGCCC', 'GCCCUUAAAACGGGGC', 'GCCCAAAAAACUGGGC', 'GGCACAGAAAAGUGCC', 'GCCCUAAUAACAGGGC', 'GGGUGCAAGACCACCC', 'GGGGCAAAGAAGCUUC', 'GGGGCAGAACAGCCCC', 'GGUGCAGAAAAGUGCC', 'GCCCCAAAAGAGGGGC']

kbest: [RNAStructure('GGGGCAAAAACGCCCC', 0.9745324428444898), RNAStructure('GGGGCAGAACAGCCCC', 0.9745537982728839), RNAStructure('GGGGCAAAGAAGCCCC', 0.9745522394520901), RNAStructure('GGGGCAGAAAAGCCCC', 0.9745541635522057), RNAStructure('GGCCCAAAAAAGGGCC', 0.974676524461917), RNAStructure('GCCCCAAUAAAGGGGC', 0.9745688330920909), RNAStructure('GCCCCAAAAAAGGGGC', 0.9745692982930593), RNAStructure('GCCCCGAAAAAGGGGC', 0.9758500630343186), RNAStructure('GGGGCAAAAAAGCCCC', 0.9745543077208556), RNAStructure('GGAGCGAAGAAGCUCC', 0.9756107700034071)]

ned_best: (0.0030683500152726695, 'GCCCCGAAAAAGGGGC')

dist_best: (-2, 'CGCGCAAAAAAGCGCG')  # -2 means the sequence adopt the target structure as its unique MFE structure

full results are saved in the file: results_20250718161920_1975123.json  # results_{timestamp}_{random_id}.json
```

### Batch mode (input a file in which each line is a structure)
``python main.py --t 1 --k 10 --object pd --path data/eterna/eterna100.txt # plus other design options`` 

### Design options
- `--online`: Enable online mode (read a single structure from stdin).
- `--t`: Sampling temperature T. Higher T increases randomness; lower T is greedier. Default: 1.
- `--k`: Frontier (priority queue) size. Larger values explore more candidates per iteration at higher cost. Default: 10.
- `--object`: Optimization objective: `pd` (probability defect) or `ned` (normalized ensemble defect). Default: `pd`.
- `--init`: Sequence initialization method: `cg` (constraint-guided/targeted) or `all` (uniform random). Default: `cg`.
- `--nosm`: Disable structured mutation (ablates structure-aware grouped changes).
- `--nomfe`: Disable “MFE as product”; only the single best-scoring sequence is kept as the designed output.
- `--step`: Maximum number of iterations. Default: 5000.
- `--repeat`: Number of repeated runs (batch mode only). Default: 1.
- `--start`: Starting index for naming output files (batch mode only). Default: 0.
- `--batch_size`: Number of structures processed concurrently in batch mode.
- `--worker_count`: Number of CPU workers for batch mode; must be ≤ `batch_size`.

### Batch mode output format
The design results will be output as a csv file with the following columns.

- ``structure``: the input puzzle (target secondary structure). 
- ``objective``: the best (lowest) objective value achived during optimization. 
- ``rna``: the sequence with the best objective during optimization. 
- ``mfe``: the MFE structure of ``rna``, given by ViennaRNA folding engine. 
- ``dist``: the distance between ``mfe`` and ``structure``. 
- ``mfe_list``: a list containing all the MFE solutions found. The order of the list is based on the time when each MFE solution is detected. 
- ``umfe_list``: a list containing all the uMFE solutions found. The order of the list is based on the time when each uMFE solution is detected.
- ``k_best``: the priority queue at the final iteration.
- ``ned_best``: the best ned value and the corresponding sequence.
- ``dist_best``: the best distance value and the corresponding sequence.
- ``log``: the objective values at each iteration.
- ``time``: the total time used to design the input puzzle ``structure``.
