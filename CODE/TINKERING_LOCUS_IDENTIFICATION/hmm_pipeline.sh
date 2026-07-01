#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1


gff="INPUT NAME OF GFF HERE"
clusters="clusters.txt" # file containing clusters of proteins in annotation proteome; contained here
contig_lengths="contig_lengths" # file containing lists of contigs in the given species (generated with esl-seqstat; format is contig_ID, space, length in nt; for simulation of randomly dispersed sequences
# this is not contained here and must be generated from each species' fastsa file as
# esl-seqstat (fasta) | awk '{print $2, $3}' > contig_lengths
# you must then manually remove the lines at the bottom with summary statistics

# make composite GFF used for analysis
grep -f conserved_extraction_keys "$gff" > "$gff"_CONSERVED # extracts gff of only conserved exons
grep -f fragment_extraction_keys "$gff" > "$gff"_FRAGMENTS # extracts gff of only fragments; for merging and then reassembly
grep -vf remove_list "$gff"_FRAGMENTS > "$gff"_FRAGMENTS_NOHIST_NOMITO # removes NUMTS and histone arrays
python merge_fragments.py "$gff"_FRAGMENTS_NOHIST_NOMITO "$clusters" # merges adjacent fragments, eg adjacent exons from the same gene, so as not to count them twice in estimating number of dispersed fragments
cat "$gff"_CONSERVED >> "$gff"_CONSERVED_PLUS_FRAGMENTS_NOHIST_NOMITO_MERGEDCLUST
cat "$gff"_CONSERVED_FRAGMENTS_NOHIST_NOMITO_MERGEDCLUST >> "$gff"_CONSERVED_PLUS_FRAGMENTS_NOHIST_NOMITO_MERGEDCLUST # reassemble gff of conserved and fragments; to identify dispersed fragments

# use the gff to find dispersed fragments and output a list of distances between them, as well as a gff with just the dispersed fragments
# requires contig lengths file and above gff
# outputs dispersed_fragment_pairwise_distances and dispersed_fragments.gff
python find_dispersed_fragments.py "$gff"_CONSERVED_PLUS_FRAGMENTS_NOHIST_NOMITO_MERGEDCLUST "$contig_lengths"

# run fit to test one vs two parameter exponential model and get best fit parameters for two parameter exponential 
# just takes above pairwise distances
python superexponential.py dispersed_fragment_pairwise_distances > dispersed_fragment_pairwise_distances_FIT
### the best fit parameters from this output file have to be read off from this output file and manually entered into the command for running the HMM below! 

# run HMM to find tinkering loci; uses gff of dispersed fragments generated above, and contig lengths file, and initialization values from the above fit
# will write run info to HMM_RUN_DETAILS and segmentation of regions to segmentation_out
python HMM.py dispersed_fragments.gff contig_lengths segmentation_out --lam_hot 2.1249e-04 --lam_cold 1.2094e-05 --max_iter 1000 --p_stay_hot 0.9997 --p_stay_cold 0.99998 > HMM_RUN_DETAILS
