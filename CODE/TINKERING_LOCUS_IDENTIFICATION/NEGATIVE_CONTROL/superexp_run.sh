#!/bin/bash
#SBATCH -J superexp_run
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-100:00
#SBATCH -c 33
python superexponential.py randomized_pairwise_distances > randomized_pairwise_distances_superexp_output
