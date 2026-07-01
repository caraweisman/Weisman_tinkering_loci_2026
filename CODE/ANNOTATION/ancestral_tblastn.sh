#!/bin/bash
#SBATCH -J tblastn
#SBATCH -n 1
#SBATCH -N 1

db="GCA_035041775.1_ASM3504177v1_genomic.fna"

makeblastdb -in "$db" -dbtype nucl

tblastn -query GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa -db "$db" -evalue 0.01 -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length frames" > "$db"_tblastn_blastout


