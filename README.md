# SAMFEO
RNA Design via Structure-Aware Multi-Frontier Ensemble Optimization
<p align="center">
<img src="figs/structure.png" width="480">
</p>

## Dependency
python3 \
ViennaRNA-2.5.1

## Environment Variable
``export VIENNAPATH=/path/to/ViennaRNA/lib/python3.9/site-packages``

## Structured Mutation
<p align="center">
<img src="figs/sm.png" width="400">
</p>
Diagrams of Structured Mutation. Paired positions are connected by blue dashed lines and each of the rounded rectangles or triangles represents a specific local structure. Diagram a,b,c show structured mutation with paired positions(shaded nucleotides pair). Diagram d,e,f show structured mutation with unpaired positions(shaded nucleotides). When a shaded position is selected for mutation, all the positions within the same local structure would be mutated simutaneously.

## Eterna100 Design
``python main.py --t 1 --k 10 --object pd --path data/eterna/eterna100.txt`` 

(Design results will be output as a csv file)
<img src="figs/inter.png" width="600">

## Base Pair Probability plot
Blue denotes pairs in the target structure and Red denotes predicted pairs not in the target structure. The darkness of the curves indicates pairing probability. The following figs show the design results for the puzzle "Runner" by NEMO and SAMFEO respectively. We can see that there are some red curves in the plot of NEMO, which means those undesired pairs are more likely to apear in the ensemble of sequence designed by NEMO.
### NEMO
<img src="figs/eterna61_nemo.png" width="800">

### SAMFEO
<img src="figs/eterna61_samfeo.png" width="800">

## 16S Design
![alt text](figs/long_design.png)
