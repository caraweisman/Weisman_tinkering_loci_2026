#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 30

query="GCA_035041775.1_ASM3504177v1_genomic.fna_esltrans"

hmmscan --acc --cpu 29 --domtblout TE_HMM_DB_vs_"$query"_hmmsearch_domtblout_ACC TE_DATABASE/TE_HMM_DB "$query" 
