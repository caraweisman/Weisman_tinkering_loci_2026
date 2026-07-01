#!/bin/bash
#SBATCH -J ccau_hmm
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-100:00
python bw_pd_2_post.py filtered_frags.gff contig_lengths segmentation_out --lam_hot 2.1249e-04 --lam_cold 1.2094e-05 --max_iter 1000 --p_stay_hot 0.9997 --p_stay_cold 0.99998 --threshold 0.5 --save_posteriors > run_output
