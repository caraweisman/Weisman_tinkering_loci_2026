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

# extract composite transcript structures from the above general output file 
# bad naming, sorry -- chimera here meant 'fragments from multiple genes on transcript'. in combination with ORF_SINGLETON, which denotes that only one of them contributes to 
# an intact ORF, this means: a composite. 
awk '$1=="CHIMERA"' allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC | grep ORF_SINGLETON > allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY
awk '$1=="CHIMERA"' allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC | grep ORF_SINGLETON > allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY

# compress into list of protein-locus pairs 
python cluster_isoforms_products.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY  > SUMMARY_OUTFILE
python cluster_isoforms_products.py allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY  >> SUMMARY_OUTFILE

# aggregate these lists over sexes
python cluster_isoforms_products_aggregate.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC_NCCHIMSONLY  >> SUMMARY_OUTFILE

