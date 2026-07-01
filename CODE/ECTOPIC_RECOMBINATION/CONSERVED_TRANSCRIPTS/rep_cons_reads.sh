#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 250000
#SBATCH -t 0-100:00
#SBATCH -c 1

femalebam="both_D_affinis_B_00_F_pacbio.fastq_refseq.sorted.bam"
malebam="both_D_affinis_B_00_M_pacbio.fastq_refseq.sorted.bam"
maingff="affinisB-merged.gff_MAINONLY"
femalefasta="both_D_affinis_B_00_F_pacbio.fasta"
malefasta="both_D_affinis_B_00_M_pacbio.fasta"

# find reads that overlap with annotated conserved genes for male/female isoseq data
bedtools intersect -abam "$femalebam" -b "$maingff" -wb -bed -split > MAIN_"$femalebam"_gffinfo
bedtools intersect -abam "$malebam" -b "$maingff" -wb -bed -split > MAIN_"$malebam"_gffinfo
# concatenate data across sexes, both bedtols output and the fasta files themselves 
cat *gffinfo >> MAIN_gffinfo
cat "$malefasta" >> all_isoseq_reads.fasta
cat "$femalefasta" >> all_isoseq_reads.fasta

# cluster reads and extract representative reads (= isoforms) for each annotated protein
# usage = find_reads_by_cons.py (bedtools-derived gff, above) (fasta) 
python -u find_reads_by_cons.py MAIN_gffinfo all_isoseq_reads.fasta

# get coordinates of features (exons) on each read so that parental exon can be excluded from analysis (to match genomic analysis) 
# use pysam to get this; desired coluns are 'read-start' and 'read_end', the coordinates of the feature on the read 
#for sex in {M,F}; do
    bam="both_D_affinis_B_00_${sex}_pacbio.fastq_refseq.sorted.bam"
    python feature_position_on_reads.py MAIN_gffinfo "$bam" "feature_positions_on_reads_${sex}.tsv"
done
 { head -n 1 feature_positions_on_reads_F.tsv; tail -n +2 feature_positions_on_reads_F.tsv; tail -n +2 feature_positions_on_reads_M.tsv; } > feature_positions_on_reads_all.tsv

# then used grep and cut to get all read ids from rep read fastas and put them all iN ALL_REP_READIDS
# then uniqued - not sure why tehy're duplicated
sort ALL_REP_READIDS | uniq > ALL_REP_READIDS_UNIQ

# pull out coordinate info only for rep reads
# faster awk command than I would usually care enough to do - intractable otherwise
LC_ALL=C awk 'NR==FNR {ids[$1]; next} $1 in ids' ALL_REP_READIDS_UNIQ feature_positions_on_reads_all.tsv > feature_positions_on_reads_REPREADS.tsv

