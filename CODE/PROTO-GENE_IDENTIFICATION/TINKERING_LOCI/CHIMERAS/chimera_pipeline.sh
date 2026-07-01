#!/bin/bash
#SBATCH -J chimera_pipeline
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-10:00
#SBATCH -c 6

orf_blastout_M="CHIMERIC_ORF_POSITION_INFO_TL_reads_M.fna_ORIENTED_m_esltrans"
orfs_M="TL_reads_M.fna_ORIENTED_m_esltrans"
reads_to_locus_M="allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads"

orf_blastout_F="CHIMERIC_ORF_POSITION_INFO_TL_reads_F.fna_ORIENTED_m_esltrans"
orfs_F="TL_reads_F.fna_ORIENTED_m_esltrans"
reads_to_locus_F="allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads"

Dath_genome="GCA_008121215.1_UCBerk_Dath_EB_1.0_genomic.fna"
Dazt_genome="GCA_005876895.1_DaztRS1_genomic.fna"
Dalg_genome="GCA_035041765.1_ASM3504176v1_genomic.fna"
Dhel_genome="GCA_963969585.1.fasta"

## male

# 1. filter out orfs with "two hits" that actually totally overlap
python remove_spurious_chimeras.py "$orf_blastout_M"

# 2. extract the peptides and blast them against outgroup genome 
cat "$orf_blastout_M"_filtered | awk '{print $1}' | uniq > CHIMERIC_ORFLIST_M
esl-sfetch --index "$orfs_M" 
esl-sfetch -f "$orfs_M" CHIMERIC_ORFLIST_M > CHIMERIC_ORFS_M.fa

makeblastdb -in "$Dath_genome" -dbtype nucl
makeblastdb -in "$Dalg_genome" -dbtype nucl
makeblastdb -in "$Dazt_genome" -dbtype nucl
makeblastdb -in "$Dhel_genome" -dbtype nucl

tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_M.fa -db "$Dath_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_M.fa_vs_"$Dath_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_M.fa -db "$Dalg_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_M.fa_vs_"$Dalg_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_M.fa -db "$Dazt_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_M.fa_vs_"$Dazt_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_M.fa -db "$Dhel_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_M.fa_vs_"$Dhel_genome"_blastout

# compare focal and target blasts to see if protein breakpoints in chimeric orf are close and in right relative orientation/order in target
python multispec_check_conserved.py "$orf_blastout_M"_filtered "$orfs_M" "$reads_to_locus_M" CHIMERIC_ORFLIST_M.fa_vs_"$Dath_genome"_blastout CHIMERIC_ORFLIST_M.fa_vs_"$Dalg_genome"_blastout CHIMERIC_ORFLIST_M.fa_vs_"$Dazt_genome"_blastout CHIMERIC_ORFLIST_M.fa_vs_"$Dhel_genome"_blastout

## female 

# 1. filter out orfs with "two hits" that actually totally overlap
python remove_spurious_chimeras.py "$orf_blastout_F"

# 2. extract the peptides and blast them against outgroup genome 
cat "$orf_blastout_F"_filtered | awk '{print $1}' | uniq > CHIMERIC_ORFLIST_F
esl-sfetch --index "$orfs_F" 
esl-sfetch -f "$orfs_F" CHIMERIC_ORFLIST_F > CHIMERIC_ORFS_F.fa

tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_F.fa -db "$Dath_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_F.fa_vs_"$Dath_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_F.fa -db "$Dalg_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_F.fa_vs_"$Dalg_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_F.fa -db "$Dazt_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_F.fa_vs_"$Dazt_genome"_blastout
tblastn -num_threads 5 -evalue 0.1 -query CHIMERIC_ORFS_F.fa -db "$Dhel_genome" -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sframe"  > CHIMERIC_ORFLIST_F.fa_vs_"$Dhel_genome"_blastout

# compare focal and target blasts to see if protein breakpoints in chimeric orf are close and in right relative orientation/order in target
python multispec_check_conserved.py "$orf_blastout_F"_filtered "$orfs_F" "$reads_to_locus_F" CHIMERIC_ORFLIST_F.fa_vs_"$Dath_genome"_blastout CHIMERIC_ORFLIST_F.fa_vs_"$Dalg_genome"_blastout CHIMERIC_ORFLIST_F.fa_vs_"$Dazt_genome"_blastout CHIMERIC_ORFLIST_F.fa_vs_"$Dhel_genome"_blastout

