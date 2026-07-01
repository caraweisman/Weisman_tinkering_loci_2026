#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

# usage: python hmmer_TEs_to_gff.py (hmmeroutfile) (orf list) (outfile name)
hmmeroutfile="TE_HMM_DB_vs_GCA_035041775.1_ASM3504177v1_genomic.fna_esltrans_hmmsearch_domtblout_ACC"
orflist="GCA_035041775.1_ASM3504177v1_genomic.fna_esltrans_ORFIDs"
outfilename="HMMER_TE_ccau.gff"

python hmmer_TEs_to_gff.py "$hmmeroutfile" "$orflist" "$outfilename"
