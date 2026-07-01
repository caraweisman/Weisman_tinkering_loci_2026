#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

# name of the gff containing the real dispersed fragments based on which you wish to perform the negative control;
# the code will use the number and length distribution of fragments in this GFF to simulate their random distribution in the genome
gff="INPUT NAME OF GFF HERE"
# similarly, contig lengths file of genome for which you wish to perform the simulation
# generated from a fasta file using the esl package in hmmer
# esl-seqstat -a (fasta) | awk '{print $2, $3}' 
# then remove the trailing summary lines
contig_lengths="contig_lengths"

# randomize the fragments
python randomize_fragments.py "$gff" "$contig_lengths"
# outputs randomized version of the gff: randomized_dispersed_fragments.gff

# compute pairwise distances between ranndomized fragments
# outputs randomized_pairwise_distances
python randomized_nn_distances.py

# run fit to test one vs two parameter exponential model and get best fit parameters for two parameter exponential 
# just takes above pairwise distances
python superexponential.py randomized_pairwise_distances > randomized_pairwise_distances_superexp_output

# run HMM to find tinkering loci; uses gff of dispersed fragments generated above, and contig lengths file, and initialization values from the above fit
### YOU HAVE TO LOOK AT THE ABOVE OUTPUT FILE, GET THE BEST FIT LAMBDA HOT AND LAMBDA COLD VALUES, AND PUT THEM BELOW
# will write run info to HMM_RUN_DETAILS and segmentation of regions to segmentation_out
python HMM.py randomized_dispersed_fragments.gff contig_lengths segmentation_out --lam_hot 2.1249e-04 --lam_cold 1.2094e-05 --max_iter 1000 --p_stay_hot 0.9997 --p_stay_cold 0.99998 > HMM_RUN_DETAILS
