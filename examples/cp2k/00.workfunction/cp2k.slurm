#!/bin/bash

#SBATCH -J cp2k
#SBATCH -o cp2k.out.%j
#SBATCH -e cp2k.err.%j
#SBATCH -p c52-small
#SBATCH -N 1
#SBATCH --ntasks-per-node=28
#SBATCH --exclusive

# add modulefiles
ulimit -s unlimited

module load miniconda/3
source activate gnn

module load intel/17.5.239 mpi/intel/2017.5.239
module load gcc/5.5.0
module load cp2k/7.1

python run.py
