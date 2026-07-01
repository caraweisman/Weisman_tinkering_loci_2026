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

# find reads that intersect with tinkering loci GFF features
bedtools intersect -abam "$femalebam" -b "$TLgff" -wb -bed -split > TL_reads_F_out_gffinfo
bedtools intersect -abam "$malebam" -b "$TLgff" -wb -bed -split > TL_reads_M_out_gffinfo

# extract read IDs for overlapping reads
cat TL_reads_F_out_gffinfo | awk '{print $4}' | sort | uniq > TL_readids_F
cat TL_reads_M_out_gffinfo | awk '{print $4}' | sort | uniq > TL_readids_M

# find reads that intersect ALL features
# also involves determining strand information correctly; TS field is often blank, since minimap can't tell transcript direction if no splice sites are present
## last three columns, in order, are: 
## bio: desired field: transcript vs reference strand. results from an operation on "fs" and "ts": if same, plus; if different, -. 
## fs and ts are derived from minimap/samtools and it's hard to extract the info; the complicated awk does that.
## fs: orientation of read relative to reference. column 6 of gff output, from sam flag 16. only differs from column 6 in cases of secondary/supp alignments (code chooses first for read). 
## ts: minimap tries to get this from splice sites. it sometimes fails. if it fails, it prints '.'
## so columns are bio, fs, ts
## bio can be . only when ts is .
## in this case, need to use polyA polarization to get ts
bedtools intersect -abam "$femalebam" -b "$fullgff" -wb -bed -split | awk -v lk=<(samtools view "$femalebam" | awk 'BEGIN{OFS="\t"} {fs=(and($2,16))?"-":"+"; ts="."; for(i=12;i<=NF;i++) if($i~/^ts:A:/) ts=substr($i,6,1); bio=(ts=="+")?fs:((ts=="-")?((fs=="+")?"-":"+"):"."); print $1,bio,fs,ts}' | sort -u -k1,1) 'BEGIN{OFS="\t"; while((getline l < lk)>0){split(l,a,"\t"); b[a[1]]=a[2]"\t"a[3]"\t"a[4]}} {print $0, (($4 in b)?b[$4]:".\t.\t.")}' > allfeature_reads_F_out_gffinfo
bedtools intersect -abam "$malebam" -b "$fullgff" -wb -bed -split | awk -v lk=<(samtools view "$malebam" | awk 'BEGIN{OFS="\t"} {fs=(and($2,16))?"-":"+"; ts="."; for(i=12;i<=NF;i++) if($i~/^ts:A:/) ts=substr($i,6,1); bio=(ts=="+")?fs:((ts=="-")?((fs=="+")?"-":"+"):"."); print $1,bio,fs,ts}' | sort -u -k1,1) 'BEGIN{OFS="\t"; while((getline l < lk)>0){split(l,a,"\t"); b[a[1]]=a[2]"\t"a[3]"\t"a[4]}} {print $0, (($4 in b)?b[$4]:".\t.\t.")}' > allfeature_reads_M_out_gffinfo

# extract reads 
esl-sfetch --index "$malefasta"
esl-sfetch --index "$femalefasta"
esl-sfetch -f "$malefasta" TL_readids_M > TL_reads_M.fna
esl-sfetch -f "$femalefasta" TL_readids_F > TL_reads_F.fna

# orient reads using polyA tail detection so that they're all the coding strand to make downstream translation simpler AND generate a table with the equivlaent of ts for when minimap cant strand properly
python reorient_reads.py TL_reads_F.fna TL_reads_F_ORIENTED.fna
python reorient_reads.py TL_reads_M.fna TL_reads_M_ORIENTED.fna

# extract lines from full gff only for reads that overlap at least one TL feature to see everything they overlap
grep -f TL_readids_F allfeature_reads_F_out_gffinfo > allfeature_reads_F_out_gffinfo_TLREADSONLY
grep -f TL_readids_M allfeature_reads_M_out_gffinfo > allfeature_reads_M_out_gffinfo_TLREADSONLY

# python script to classify features 
python dev_find_chimeric_reads.py allfeature_reads_M_out_gffinfo_TLREADSONLY TL_reads_M.fna_TS_STRANDS 
python dev_find_chimeric_reads.py allfeature_reads_F_out_gffinfo_TLREADSONLY TL_reads_F.fna_TS_STRANDS

# esl-translate reads to find orfs; since reoriented, only translate plus strand
esl-translate -m --watson TL_reads_M_ORIENTED.fna > TL_reads_M.fna_ORIENTED_m_esltrans
esl-translate -m --watson TL_reads_F_ORIENTED.fna > TL_reads_F.fna_ORIENTED_m_esltrans

# search reads to see if they actually encode an intact orf (esl-trans -m) matching one or more of the nt components on the read
# this takes about an hour (for both M and F)!
python -u 4protein_search.py allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads TL_reads_M.fna_ORIENTED_m_esltrans "$refproteome"
python -u 4protein_search.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads TL_reads_F.fna_ORIENTED_m_esltrans "$refproteome"

# extract rows with proteins that are only present in affinis based on conservation analysis
python 2coding_conservation.py allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS "$consfile" PROT_LENGTH_INFO_TL_reads_M.fna_ORIENTED_m_esltrans
python 2coding_conservation.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS "$consfile" PROT_LENGTH_INFO_TL_reads_F.fna_ORIENTED_m_esltrans

# compress into list of protein-locus pairs
python cluster_isoforms_products.py allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC  > SUMMARY_OUTFILE
python cluster_isoforms_products.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC  >> SUMMARY_OUTFILE

# aggregate compression across sexes
python cluster_isoforms_products_aggregate.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS >> SUMMARY_OUTFILE
python cluster_isoforms_products_aggregate.py allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC allfeature_reads_M_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC  >> SUMMARY_OUTFILE

## SUMMARY OUTFILE CONTAINS SUMMARY STATISTICS ON NUMBER OF PROTEINS AND LOCI PER SEX AND FOR BOTH SEXES


