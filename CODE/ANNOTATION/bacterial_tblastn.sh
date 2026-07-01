#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 3

target="GCA_035041775.1_ASM3504177v1_genomic.fna"

makeblastdb -in "$target" -dbtype nucl

for file in *.faa; do tblastn -query "$file" -db "$target" -evalue 0.01 -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length frames" -num_threads 2 >  "$file"_vs_"$target"_annotation_tblastn_blastout; done

