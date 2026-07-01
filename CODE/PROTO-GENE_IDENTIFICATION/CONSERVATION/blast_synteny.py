# updated 5/20 to make compatible with new 5/20 synteny0 script which potentially outputs L and R flank synteny coordinates per fragment, named (eslid)-R/L

import sys
import os
import re

focal_species_fasta_file = sys.argv[1]
target_species_fasta_file = sys.argv[2]

os.system('esl-sfetch --index ' + focal_species_fasta_file)
os.system('esl-sfetch --index ' + target_species_fasta_file)

focal_ids = set()
for line in open(focal_species_fasta_file):
	if line.startswith('>'):
		focal_id = re.split(' ', line)[0][1:]
		focal_ids.add(focal_id)

target_ids = set()
for line in open(target_species_fasta_file):
	if line.startswith('>'):
		target_id = re.split(' ', line)[0][1:]
		target_ids.add(target_id)

common_ids = []
for focal_id in focal_ids: 
	for target_id in target_ids:
		if target_id == focal_id or target_id.startswith(focal_id + '-'): 
			common_ids.append([focal_id, target_id])


# clear from previous runs
os.system('rm ' + focal_species_fasta_file + '_' + target_species_fasta_file + '_blastout')

for seqids in common_ids:
	focalid = seqids[0]
	targetid = seqids[1]
	eslcmd1 = 'esl-sfetch ' + focal_species_fasta_file + ' ' + focalid + ' > query.fna'
	eslcmd2 = 'esl-sfetch ' + target_species_fasta_file + ' ' + targetid + ' > target.fna'
	print(eslcmd1)
	print(eslcmd2)
	os.system(eslcmd1)
	os.system(eslcmd2)
	dbcmd = 'makeblastdb -in target.fna -dbtype nucl'
	os.system(dbcmd)
	print(dbcmd)
	blastcmd = 'blastn -task dc-megablast -evalue 0.5 -num_threads 5 -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue score pident length sstrand qcovs" -query query.fna -db target.fna >> ' + focal_species_fasta_file + '_vs_' + target_species_fasta_file + '_blastout'
	print(blastcmd)
	os.system(blastcmd) 

