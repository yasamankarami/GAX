# GAX
Genetic Algorithm for minimum ensemble fitting of SAXS data


## Installation using conda:

conda create --yes -n GAX -c conda-forge python=3.8 numpy scipy MDAnalysis MDAnalysisTests

source activate GAX

conda install imp


## Running GAX:

python GA_wrapper.py --addr PATH_TO_THE_INPUT_FILES --pdb TOPOLOGY --trj TRAJECTORY --saxs EXPERIMENTAL_SAXS_DATA -s 5 10 -r 2
