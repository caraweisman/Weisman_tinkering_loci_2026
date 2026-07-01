#!/bin/bash
#SBATCH -J mmseqs_cluster
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-100:00
#SBATCH --cpus-per-task 31

# Step 1: make a DB from your reads
mmseqs createdb F_unannotated_reads.fna F_unannotated_reads.fnaDB

# Step 2: cluster the DB. this clusters sequences if they are 100% identical over the max of the two sequences' lengths, ie exact deduplication. 
mmseqs cluster F_unannotated_reads.fnaDB clusterDB tmp --min-seq-id 1 -c 1 --cov-mode 0 --threads 31

# Step 3: extract representatives
mmseqs result2repseq F_unannotated_reads.fnaDB clusterDB repDB

# Step 4: convert back to FASTA
mmseqs convert2fasta repDB F_unannot_transcripts_uniq.fa

# cluster membership table
mmseqs createtsv F_unannotated_reads.fnaDB F_unannotated_reads.fnaDB clusterDB cluster_membership.tsv

# counts per representative
cut -f1 cluster_membership.tsv | sort | uniq -c | sort -nr > cluster_sizes.txt


