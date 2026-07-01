#!/bin/bash
#SBATCH -J DALG_FRAGCONS
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p eddy
#SBATCH --mem 150000
#SBATCH --exclude=holygpu7c0920,holygpu8a2650[4-6]
#SBATCH -t 0-10:00
#SBATCH -c 31


##### INPUT FILES

# string for species identifier; can be replaced with anything
species="Dalg"
# focal species GFF; switch for tinkering loci vs outside tinkering loci. Both provided in the GFF directory.
focal_gff="TLfrags_PLUSMAIN.gff"
# target species GFF; combination of the conserved genes plus fragments. Not provided here again (too much space for all 5 outgroups), but can be produced by
# following the instructions in and using the original files in the ANNOTATIONS supplementary files. 
target_gff="dalg_merged.gff_MAINPLUSFRAG"
# contig lengths file for affinis; provided here
focal_contig_lengths="contig_lengths"

# affinis genome; download from NCBI (accession in Table S6)
focal_genome="GCA_035045985.1_ASM3504598v1_genomic.fna"
# target species genome; download from NCBI (Accession in Table S6)
target_genome="GCA_035041765.1_ASM3504176v1_genomic.fna"




## run step 0 of synteny search
# this 1. identifies the flanking conserved genes on either side of each fragment in affinis. 
# this is either a) the region between the two flanking conserved genes
# b) if the fragment is next to a contig end, the region between the nearest conserved exon and that end
# or c) if the fragment is on its own contig with no conserved genes, the entirety of that contig
# then 2. it prints a list of these coordinates (in a file that can be used as input with esl-sfetch -Cf) for ALL genes to focal_fragment_synteny_eslcoords.
# then 3. for only case a) above, where there are TWO marker genes so there is no ambiguity, 
# it determines whether those two genes are also direct neighbors in the outgroup. 
# if they are, it looks for a fragment with the same annotation in the outgroup. for each fragment in the focal species, 
# it prints True/False for this search to 0_syntenyfound_list. (if no properly flanked syntenic region found in focal species, this is automatically False.) 
# if there is an intact corresponding syntenic region in the outgroup with no identical fragment annotation, it prints the genomic coordinates of the outgroup region
# to target_fragment_synteny_eslcoords. 
python -u synteny0.py "$focal_gff" "$target_gff" "$focal_contig_lengths"
# outputs 1. 0_syntenyfound_list, the True/False assessment above
# 2. focal_fragment_synteny_eslcoords, list of coordinates in focal species (all fragments)
# 3. target_fragment_synteny_eslcoords, list of coordinates in target species (fragments with intact 2 flanking conserved genes)
# 4. synteny_broken_list, list of fragments where there are 2 flanking conserved genes in focal species that are not adjacent in target (for interest

## extract esl coordinates of all focal fragments that did not have a syntenic fragment identified in the above approach
# (output them all in the first place just because it might be useful; not necessarily going to use them later)
# uniq is because some fragments are twice in gff. this is annoying and I don't know why - I think because of variant isoforms that are smaller than my deoverlap threshold
grep False 0_syntenyfound_list | awk '{print $1}' | grep -Ff - focal_fragment_synteny_eslcoords | sort | uniq > 0unfound_focal_fragment_synteny_eslcoords

# another dedup
cat target_fragment_synteny_eslcoords | sort | uniq >  target_fragment_synteny_eslcoords2
mv target_fragment_synteny_eslcoords2 target_fragment_synteny_eslcoords

## extract the focal and target syntenic regions as fasta sequences using esl 
esl-sfetch --index "$focal_genome"
esl-sfetch --index "$target_genome"
esl-sfetch -Cf "$focal_genome" 0unfound_focal_fragment_synteny_eslcoords > 0unfound_Daff_synteny_regions.fna
esl-sfetch -Cf "$target_genome" target_fragment_synteny_eslcoords > "$species"_synteny_regions.fna

## run step 1 of synteny search 
# run blast to compare focal syntenic regions with syntenic gene-identified syntenic regions in outgroup
# use blast for this because there is a strong a priori hypothesis about where they should be
# therefore can just look for whether the fragment is there at all within that region; no reconstruction around the highly perturbable/fragmented blast results required
# include focal flanking sequence to sanity check to see whether it actually looks at all orthologous (vs something weird happening, like a huge insertion from another region that should then be found
# and checked) 
# for runtime, the script submits jobs one at a time, comparing only the same syntenic region in each species
# but outputs to a single file
#! NB this script includes parsing of sequence names that looks for identical names in the esl-sfetch-derived fasta files by cutting off the contig information from
# the source fastas, which is separated from the new eslid (in the esl sfetch file, first column) by a space. not sure if this will always work. 
python -u blast_synteny.py 0unfound_Daff_synteny_regions.fna "$species"_synteny_regions.fna
# outputs focal_species_fasta_file + '_' + target_species_fasta_file + '_blastout' : blastout file for the search

# use blast results to assess a) whether the fragment is in the syntenic region
# this also outputs b) whether there was a syntenic region with flanking conserved genes identified in the target species at all
# (this is contained in the information above, but is a useful summary, as if the answer to b) is no, the answer to a) must be no) 
# and c) if there was, whether that region covers at least 25% of the focal species' region  
python check_synteny_blast.py 0unfound_Daff_synteny_regions.fna_vs_"$species"_synteny_regions.fna_blastout 0unfound_focal_fragment_synteny_eslcoords "$species"
# outputs synteny1_results_' + species, a file with the following columns: 
# 1. frag ID; 2. whether the syntenic region was identified in the target species (b); 3. whether a syntenic region is >25% similar to the region in the focal species (c); 
# 4. whether the fragment is in the syntenic region (a)
 
# for fragments that did not have an intact syntenic region in the outgroup (c)
# or with an intact syntenc region that does not appear to really match the focal species (b)
# do a global alignment with minimap
# first filter the results from the blast synteny check to outputs IDs of a) regions with no syntenic region identified in outgroup
# and b) regions that neither have the fragment found in that putative region AND that region is such that it shares little sequence 
# with the focal syntenic region, raising suspicion that the real locus is elsewhere
# columns in output are
#Frag_ID         Syntenic_region_present_outgroup        Syntenic_region_matched_outgroup        Synteny_true 
cat synteny1_results_"$species" | awk '($2=="False") || ($3=="False" && $4=="False")' | awk '{print $1}' | sort | uniq > "$species"_syn2_globalsearch_ids

# extract these fragment names and get sequences of their syntenic regions in the focal species
esl-sfetch -f 0unfound_Daff_synteny_regions.fna "$species"_syn2_globalsearch_ids > "$species"_syn2_globalsearch_seqs.fna

## Run step 2 of synteny search
# run minimap to see if there's the fragment with at least 3 kb on one side anywhere else in the genome
minimap2 -g 25 -r 12 -m 0 -n 20 -p 0.1 -t 30 "$target_genome" "$species"_syn2_globalsearch_seqs.fna > "$species"_syn2_globalsearch_seqs.fna_vs_"$target_genome"_mm_globalsearch.paf

# then run a script that determines whether the minimap search found the orthologous fragment
# criterion is that it plus some amount of flanking sequence (currently 5 kb) are found contiguously
# minimap is good for this because it bridges gaps that blast would otherwise break up, making this easy 
python -u check_globalsynteny_minimap.py "$species"_syn2_globalsearch_seqs.fna_vs_"$target_genome"_mm_globalsearch.paf 0unfound_focal_fragment_synteny_eslcoords "$species"_syn2_globalsearch_ids "$species"
# outputs 2globalsynteny_results_' + species: 
# column 1: frag iD; column 2: whether the fragment was found by the above criteria. 

## concatenate results from all three steps
grep True 0_syntenyfound_list  >> final_syntenyfound_list 
awk -v OFS='\t' '$4=="True" {print $1, $4}' synteny1_results_"$species" >> final_syntenyfound_list
grep True 2globalsynteny_results_"$species" >> final_syntenyfound_list
sort final_syntenyfound_list | uniq > final_syntenyfound_list2
mv final_syntenyfound_list2 final_syntenyfound_list

## put into a table that lists explicitly true/false for synteny in this species
python compile_final_list.py "$species"
cat synteny_results_table_"$species" | uniq > synteny_results_table_"$species"2
mv synteny_results_table_"$species"2 synteny_results_table_"$species"
