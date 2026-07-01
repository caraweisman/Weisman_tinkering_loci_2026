#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --cpus-per-task 6

file="dfam-fasta-download.fasta"
target="GCA_035041775.1_ASM3504177v1_genomic.fna"

makeblastdb -in "$target" -dbtype nucl

blastn -task dc-megablast -query "$file" -db "$target" -evalue 0.01 -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length frames" -num_threads 5 >  "$file"_vs_"$target"_annotation_blastn_blastout

