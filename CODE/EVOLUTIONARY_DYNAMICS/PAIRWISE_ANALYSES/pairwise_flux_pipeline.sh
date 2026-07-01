#!/bin/bash
#SBATCH -J TL_flux
#SBATCH -n 1
#SBATCH -N 1

# fasta file with sequences of whole tinkering loci from focal species (affinis); in this directory
focal_TL_sequences="TLs.fna"
# a list of coordinates of each tinkering locus; used to generate the fasta file above using esl-sfetch and the genome
TL_list="TL_sfetch" # use sfetch file as list because easier than fasta and you have it anyway to make the fasta ; in this directory
# target species genome (other affinis group member, e.g. athabasca); available from NCBI
target_genome="GCA_035041765.1_ASM3504176v1_genomic.fna"
# focal species genome (here affinis); available from NCBI
focal_genome="GCA_035045985.1_ASM3504598v1_genomic.fna"
# manual label for species name; used in output file names
species="Dalg"
# annotation in GFF of focal species (affinis) tinkering locus framgents; in this directory
focal_TL_gff="TL_frags.gff"
# annotation in GFF format from target species (eg athabasca) containing only fragments; can be generated from annotation files in above SUPP_INFO directory using instructions there
target_fragment_gff="dalg_merged.gff_FRAGMENTSONLY"

# minimap to search for similar sequences to the focal TL in the outgroup
# with these parameters, should handle indels of ~200kb and call it one alignment
# so only big movements of regions will break it up in the outgroup
minimap2 -x asm20 -k 14 -w 5 -m 20 -n 2 -t 63 "$target_genome" "$focal_TL_sequences" > "${target_genome}_vs_${focal_TL_sequences}_mm.paf"; \

# use a python script to find regions that are unambiguously orthologous to the TL in outgroups
# that is, if there are mulitple regions that minimap aligns to the focal region, it does not assert an orthologous locus
# takes only subregions of each TL that have a unique mapping in the outgroup
python reconstruct_TLs.py "${target_genome}_vs_${focal_TL_sequences}_mm.paf" "$TL_list" "$species"

# do sfetch with both target and focal corresponding regions
esl-sfetch --index "$focal_genome"
esl-sfetch --index "$target_genome"
esl-sfetch -Cf "$focal_genome" TL_reconstruction_esl_focal > Daff_reconstruction_seqs.fna
esl-sfetch -Cf "$target_genome" TL_reconstruction_esl_"$species" > "$species"_reconstruction_seqs.fna

# blast them
makeblastdb -in "$species"_reconstruction_seqs.fna -dbtype nucl
blastn -num_threads 5 -task dc-megablast -query Daff_reconstruction_seqs.fna -db "$species"_reconstruction_seqs.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sstrand qcovs" > Daff_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout

# filter for self hits (sequences from query and target have same name) 
cat Daff_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout | awk '$1==$5' > Daff_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout_selfonly

# self blast on focal to detect in-duplicates
makeblastdb -in Daff_reconstruction_seqs.fna -dbtype nucl
blastn -num_threads 5 -task dc-megablast -query Daff_reconstruction_seqs.fna -db Daff_reconstruction_seqs.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sstrand qcovs" > Daff_reconstruction_seqs.fna_vs_Daff_reconstruction_seqs.fna_blastout

# filter for self hits on focal AND filter OUT the diagonal self hits
cat Daff_reconstruction_seqs.fna_vs_Daff_reconstruction_seqs.fna_blastout | awk '$1==$5 && !($3==$7 && $4==$8)' > Daff_reconstruction_seqs.fna_vs_Daff_reconstruction_seqs.fna_blastout_selfonly_nodiag

# self blast on target to detect in-duplicates
# already made blastdb 
blastn -num_threads 5 -task dc-megablast -query "$species"_reconstruction_seqs.fna -db "$species"_reconstruction_seqs.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length sstrand qcovs" > "$species"_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout

# filter for self hits on target AND filter OUT the diagonal self hit
cat "$species"_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout | awk '$1==$5 && !($3==$7 && $4==$8)' > "$species"_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout_selfonly_nodiag

# run script that finds both brand new fragments in a region and rearrangements of existing fragments
# rearranemnets includes order changes and duplications
python find_new_fragments.py Daff_reconstruction_seqs.fna_vs_"$species"_reconstruction_seqs.fna_blastout_selfonly TL_reconstruction_esl_focal TL_reconstruction_esl_"$species" "$focal_TL_gff" "$target_fragment_gff" "$species"

# output coverage of focal region in esl coordinates for phylogenetic analysis 
awk '{k=$1; sub(/-[0-9]+$/,"",k); if(!(k in s)){o[++n]=k; split(k,a,":"); split(a[2],b,"-"); d[k]=b[2]-b[1]; s[k]=1} num[k]+=$3-$2} END{for(i=1;i<=n;i++) printf "%s\t%.6f\n", o[i], num[o[i]]/d[o[i]]}' TL_reconstruction_esl_focal > focal_TL_coverage_"$species"

