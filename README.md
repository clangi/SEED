This repository is no longer mantained, as the project has been moved to the new homepage at: [https://gitlab.com/CaflischLab/SEED](https://gitlab.com/CaflischLab/SEED)

# SEED (Solvation Energy for Exhaustive Docking) version 4.0.0 

# ! DEPRECATED PLEASE VISIT [https://gitlab.com/CaflischLab/SEED](https://gitlab.com/CaflischLab/SEED) for the up-to-date version!

SEED is a program for fragment docking with force-field based evaluation of binding energy.

Please refer to the documentation (seed_4.0.0_doc.pdf) to understand how to set up a calculation in SEED.

### Installation ###
cd to directory src and make SEED with the following command (you may have to modify the Makefile.local):
```sh
make seed
```
The binary is compiled into the bin directory.
### Run SEED ###
Run SEED with the following command:
```sh
./seed_4.0.0 seed.inp >& log
```
You can find an example seed.inp and seed.par, along with the results of simple study cases,
in the test_cases folder.

In the scripts directory you can find separate_poses.py, that can be used to extract selected poses from the concatenated
mol2 output files. For using it type:
```sh
python separate_poses.py -h
```
### Citation ###
Kindly reference the original paper if you use SEED:
 * N. Majeux, M. Scarsi, J. Apostolakis, C. Ehrhardt, and A. Caflisch. Exhaustive docking of
molecular fragments on protein binding sites with electrostatic solvation.
Proteins: Structure, Function and Genetics, 37:88-105, 1999.

The description of the fast energy evaluation is in the second SEED paper:
 * N. Majeux, M. Scarsi, and A. Caflisch. Efficient electrostatic solvation model for protein-
fragment docking.
Proteins: Structure, Function and Genetics, 42:256-268, 2001.
