#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

main="GCA_035041775.1_ASM3504177v1_genomic.fna_tblastn_blastout_main_hits_HR"
extra="GCA_035041775.1_ASM3504177v1_genomic.fna_tblastn_blastout_extra_hits_pars_HR"
outfilename="ccau.gff"

hmmgff="HMMER_TE_ccau.gff"
mergedoutfilename="ccau_merged.gff"

python annots_to_gff.py "$main" "$extra" "$outfilename"
python merge_hmmerte_main_cgpt.py "$outfilename" "$hmmgff" "$mergedoutfilename"
