# GAX
Genetic Algorithm for minimum ensemble fitting of SAXS data


## Installation using conda:

conda create --yes -n GAX -c conda-forge python=3.8 numpy scipy MDAnalysis MDAnalysisTests

source activate GAX

conda install imp


## Running GAX:

If you have the back-calculated SAXS profile from the MD simulation:

usage: GA_wrapper.py [-h] --pdb PDB --saxs SAXS --mode MODE [--trj TRJ] [--calc_saxs CALC_SAXS]
             [-s S] [-r R] [-nb_gen NB_GEN] [-nb_ens NB_ENS] [-window WINDOW]

Extract ensemble of conformations matching best with the experimental SAXS data

arguments:
  -h, --help      show this help message and exit
  --pdb PDB             Protein structure file
  --saxs SAXS           Experimental SAXS data
  --mode MODE           Do you have the back-calculated profile
  --trj TRJ             Trajectory file
  --calc_saxs CALC_SAXS             Back-calculated SAXS profile

options:
  -s S                  Ensemble size
  -r R                  Number of repeats
  -nb_gen    NB_GEN                  Number of generations
  -nb_ens    NB_ENS                  Number of ensembles
  -window    WINDOW                  Window size
