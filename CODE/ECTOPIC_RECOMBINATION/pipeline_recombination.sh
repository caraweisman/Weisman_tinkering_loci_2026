#!/bin/bash
#SBATCH -J recombination_pipeline
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-100:00
#SBATCH -c 1

fraggff="TL_frags.gff" # included here
maingff="affinisB-merged.gff_MAINONLY" # included here and also derivable from files/instructions in ANNOTATION directory
contiglengths="contig_lengths" # included here
flank="2500" # adjustable parameter; size of flanking sequence
genome="GCA_035045985.1_ASM3504598v1_genomic.fna" # must be downloaded from NCBI (too large) 

#1. takes each TL fragment plus FLANK nt on either side  and blasts it against each exon in its parent gene (based on annotation) plus FLANK nt on either side
# also does a control: for each exon in the parent gene, picks an exon from a random gene and does th same thing
# outputs results to a compiled blast output file for both experimetnal and control
# also outputs coordinate conversion files: what the absolute coordinates of each segment are in the genome
# used later in combination with gff to check for TEs etc
python control_standalone_check_flank.py "$fraggff" "$maingff" "$contiglengths" "$flank" "$genome"

# 2. analyzes the blast outputs from above
# looks for hits from each fragment-parental gene pair that pass certain criteria
# hit has to be outside of the coding region (in first/last FLANK nt)
# and has to be longer than some minimum threshold length (currently 100) - idea is to find regions that are likely due to or sufficient to cause ER
# outputs a file that is a list of these successful filtered hits : recombination_hits, used in a lot of stuff
# also outputs this for the control: sites that are just as similar between frag/source putative recomb hits and between frag/random exon, which presumably did not ***** 
python regionspec_assess_blastout.py "$fraggff" "$fraggff"_source_blastout "$flank" "$fraggff"_fragment_blastout_coords "$fraggff"_source_blastout_coords "$fraggff"_CTRL_source_blastout "$fraggff"_CTRL_source_blastout_coords

# to assess relative probability of RT vs EC
# requiers preclustering of rep reads for each protein isoform, in the directory in teh code
# does blast between transcripts and frag+flanking region for regions of similarity outside of fragment
# very similar to check flank; outputs a blast file between each fragment and all representative transcripts that have been assembled for each isoform in a separate directory
# output file is "$fraggff"_TX_blastout
python check_transcript.py "$fraggff" "$flank" "$contiglengths" "$genome"

# file to output file like recombination_hits, TX_hits, showing hits from each fragment's FLANK flanking nt to noncoding regions in its transcripts; 
# blastout from above frag vs tx is input, hard-coded
# requires a minimum size of match (currently 100 nt like above) - no pident because they're old (same idea as 100 nt which is probbly too small)
# outpt file is TX_hits
# can be compared to recombination_hits to see how many fragments are supported by recombination vs RT; can do this heuristically by coutning number in output or with the below script
#python simple_tx_assess.py
python fulltx_simple_tx_assess.py

# more formally compare results from tx assessment and recombination assessment to see which is favored / if all recomb hits are explained by tx
# under the hypothesis that the shared regions are due to noncoding exons
# outfile is Compiled_recomb_TX_info; lists at the individual fragment level: 
# is there a recomb site in tx? in source? in both? in neither? 
# so you can answer - for how many fragments is ER a strictly better explanation
python compare_tx_recomb_hits.py

