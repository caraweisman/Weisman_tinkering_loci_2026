#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 100000
#SBATCH -t 0-1:30
#SBATCH -c 1

# gff of tinkering locus fragments (included here)
fraggff="TL_frags.gff"
# annotations of conserved genes (included here)
maingff="affinisB-merged.gff_MAINONLY"
# contig lengths from fasta file for genome (included here) 
contiglengths="contig_lengths"
# size of flanking sequence; 2500 used in manuscript
flank="2500"

python -u TE_enrichment.py TRANSPOSON_ANNOTATIONS/affinis_RNA_transposons_PFAM_annotations.gff DFAM "$contiglengths" "$flank" "$maingff" "$fraggff" > OUTSTATS_DFAM_ANNOTATION
python -u TE_enrichment.py TRANSPOSON_ANNOTATIONS/affinis_RNA_transposons_PFAM_annotations.gff  RNAPROT "$contiglengths" "$flank"  "$maingff" "$fraggff" > OUTSTATS_RNA_ANNOTATION

