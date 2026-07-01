#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 100000

# annotation for species in GFF format (in SUPP FILES/ANNOTATIONS directory) 
fullgff="D_affinis.gff"
# annotation just for tinkering locus fragments (in this directory)
# to analyze tinkering locus fragments, use TL_frags.gff; to analyze fragments in the rest of the genome, use outsideTL_frags.gff.
TLgff="TL_frags.gff"
# bam files from Iso-Seq aligned to reference genome, one for each sex (Iso-Seq and genome available on NCBI)
femalebam="D_affinis_B_00_F.bam"
malebam="D_affinis_B_00_M.bam"
# Iso-Seq reads in FASTA format 
femalefasta="D_affinis_B_00_F.fasta"
malefasta="D_affinis_B_00_M.fasta"
# Drosophila melanogaster reference proteome (available from NCBI)
refproteome="GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa"
# file of conservation status for fragments (in this directory) 
consfile="conservation_results"

# extract singleton loci from the above general output file 
awk '$1=="SINGLETON"' allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC > allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY
awk '$1=="SINGLETON"' allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC > allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY

# compress into list of protein-locus pairs 
# do this for both daff specific coding orfs and for all coding orfs
python cluster_isoforms_products.py allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY  > SUMMARY_OUTFILE
python cluster_isoforms_products.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY  >> SUMMARY_OUTFILE

# aggregate these lists over sexes
python cluster_isoforms_products_aggregate.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY_iso_collapsed allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_SINGLETONSONLY_iso_collapsed >> SUMMARY_OUTFILE

