## STAT525-SISPER

#### Developers:
Zhikun Cai, Chenchao Shou, Guanfeng Gao

#### Problem:
Searching for the lowest-energy folded state of protein in a 2D hydrophobic-hydrophilic (HP) lattice model

#### Algorithm:
Sequential Importance Sampling with Pilot-Exploration Resampling (SISPER)

#### Reference:
J. L. Zhang and J. S. Liu, A new sequential importance sampling method and its application to the two-dimensional hydrophobic-hydrophilic model, *Journal of Chemical Physics*, 117, 3492 (2002)

#### Code Usage:
(1) Use sisper.py to run the simulation with a prepared sequence file as input (tau acts like temperature):
```bash
    $ cd src
    $ python sisper.py sequence_file_name conformation_file_name tau
```
(2) Use plot_conformations.py to plot the protein conformation found and save the figures:
```bash
    $ python plot_conformations.py sequence_file_name conformation_file_name fig_file_keywords num_of_figs
```
    
#### Example:
See the test folder
```bash
    $ cd test
    $ python ../src/sisper.py sequence.txt conformations_tau0.5 0.5
    $ python ../src/plot_conformations.py sequence.txt conformations_tau0.5 plots_tau0.5 8
```
    
