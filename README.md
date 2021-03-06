# GAX: Genetic Algorithm for minimum ensemble fitting of SAXS data

GAX is a method based on genetic algorithm for extracting a minimum ensemble of conformations from the pool of conformations generated by molecular dynamics simulations that fit SAXS data. GAX is implemented as a webserver and an open-source python-based software.

## Installation using conda:

conda create --yes -n GAX -c conda-forge python=3.8 numpy scipy MDAnalysis MDAnalysisTests

source activate GAX

conda install imp


## Running GAX:
usage: GA_wrapper.py [-h] --pdb PDB --saxs SAXS --mode MODE [--trj TRJ] [--calc_saxs CALC_SAXS]
             [-s S] [-r R] [-nb_gen NB_GEN] [-nb_ens NB_ENS] [-window WINDOW]

GAX has two modes:

- Using back-calculated SAXS profile from a molecular dynamics simulation (--mode yes). This file has to contain a SAXS profile for every conformation.

- Using a trajectory file from a molecular dynamics simulation (--mode no). The acceptable trajecotry formats are those readable by "MDAnalysis".



# Arguments: 
  -h, --help             show this help message and exit
  
  --pdb PDB              Protein structure file
  
  --saxs SAXS            Experimental SAXS data
  
  --mode yes/no          Do you have the back-calculated profile
  
  --trj TRJ              Trajectory file
  
  --calc_saxs CALC_SAXS  Back-calculated SAXS profile

# Advanced options:

  -s S                  Ensemble size
  
  -r R                  Number of repeats
  
  -nb_gen    NB_GEN     Number of generations
  
  -nb_ens    NB_ENS     Number of ensembles
  
  -window    WINDOW     Window size
