#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

# compile coverage statistics of each affinis tinkering locus in all outgroups into a single file for ease of reading and use; filenames are hard-coded
# output of this file is included here for ease of reference/as supplemental material: compiled_coverage
python compile_coverage_stats.py
# use this to make a list of regions that are excluded because their coverage is <50% for all species
awk 'NR==1 || ($2<0.5 && $3<0.5 && $4<0.5 && $5<0.5)' compiled_coverage | awk '{print $1}' > insufficient_coverage

# take in all pairwise analysis files and calculate fragment gain flux for each tinkering locus across whole phylogeny
python calculate_flux_fragment_gain.py

# extract tinkering loci that have sufficient coverage
grep -vf insufficient_coverage Total_gainedfrags_flux_numbers > Total_gainedfrags_flux_numbers_sufficient_coverage
grep -vf insufficient_coverage Gainedfrags_flux_IDs > Gainedfrags_flux_IDs_sufficient_coverage

# take in all pairwise analysis files and calculate fragment gain flux for each tinkering locus across whole phylogeny
python calculate_flux_fragment_dup.py

grep -vf insufficient_coverage Total_duplicatedfrags_flux_numbers > Total_duplicatedfrags_flux_numbers_sufficient_coverage
grep -vf insufficient_coverage Dupfrags_flux_IDs > Dupfrags_flux_IDs_sufficient_coverage

