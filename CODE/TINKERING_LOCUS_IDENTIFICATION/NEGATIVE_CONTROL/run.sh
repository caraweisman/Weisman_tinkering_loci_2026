#!/bin/bash
#SBATCH -J hmm
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-100:00
python HMM.py randomized_dispersed_fragments.gff contig_lengths segmentation_out --lam_hot 2.3101e-05 --lam_cold 2.3101e-05 --p_stay_hot 0.9997 --p_stay_cold 0.99998 --max_iter 1000 --threshold 0.5 --save_posteriors  > RUN_INFO
