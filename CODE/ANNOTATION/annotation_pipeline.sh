#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

main_blast="GCA_035041775.1_ASM3504177v1_genomic.fna_tblastn_blastout"
bacterial_blast="Acetobcter_Lactiplanti_Wolbachia_vs_GCA_035041775.1_ASM3504177v1_genomic.fna_annotation_tblastn_blastout"
TE_blast="dfam-fasta-download.fasta_vs_GCA_035041775.1_ASM3504177v1_genomic.fna_annotation_blastn_blastout"
protein_list="Dmel_accession_info"

# find main hits; output file is input file name plus '_main_hits'
python -u find_main_hits.py "$main_blast" "$protein_list"
# find extra hits; syntax is Usage: find_extra_hits.py <full_blast> <bacterial_blast> <TE blast> <main_hits>
python -u find_extra_hits.py "$main_blast" "$bacterial_blast" "$TE_blast" "${main_blast}_main_hits"
# merge clusters; residual 
python merge_clusters.py "$main_blast"_extra_hits "$main_blast" clusters.txt 
# make new ("pars" - parsimonious) extra hits file and main hits file "human readable" by substituting acccessions for full descriptions
python human_readable_annotation.py "$main_blast"_main_hits All_fullinfo
python human_readable_annotation.py "$main_blast"_extra_hits_pars All_fullinfo


