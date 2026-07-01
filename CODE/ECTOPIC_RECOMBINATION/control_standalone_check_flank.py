import numpy as np 
import sys
import os
import re 
import subprocess
import random

fragff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)
maingff = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t', comments=None)
contig_lengths = np.genfromtxt(sys.argv[3], dtype=str)
flank = int(sys.argv[4])
genome = sys.argv[5]

os.system('esl-sfetch --index ' + genome)

def extract_acc(desc): 
	acc = re.split('-', desc)[2]
	return acc

def extract_pos(desc):
	match = re.search(r'\((\d+-\d+)/', desc)
	result = match.group(1)  # '1-652'
	return result 

contig_lengths_dict = {}
for contig in contig_lengths: 
    contigname = contig[0]
    contiglen = int(contig[1])
    contig_lengths_dict[contigname] = contiglen

main_pos_dict = {} ## format: key = accession, value = [[start, stop, contig]] for each exon
main_pos_list = [] # for random sampling of individual exons; list instead of dict because don't want to bias by number of exons
for row in maingff: 
	desc = row[8]
	acc = extract_acc(desc) 
	start = int(row[3])
	stop = int(row[4]) 
	contig = row[0]
	pos = extract_pos(desc)
	main_pos_list.append([start, stop, contig, acc])
	if acc not in main_pos_dict: 
		main_pos_dict[acc] = [[start, stop, contig, pos]]
	elif acc in main_pos_dict: 
		main_pos_dict[acc].append([start, stop, contig, pos])


os.system('rm ' + sys.argv[1] + '_source_blastout') # clear old blast outputs
os.system('rm ' + sys.argv[1] + '_CTRL_source_blastout') # clear old blast outputs
query_outfile = open(sys.argv[1] + '_fragment_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later  
target_outfile = open(sys.argv[1] + '_source_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later  
#ctrl_outfile =  open(sys.argv[1] + '_CTRL_source_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later 
#not actually necessary because coordinates are included in target id; see next script

for i in range(0, len(fragff)):
	numexons = 0
	row = fragff[i] 
	contig = row[0]
	fragstart = int(row[3])
	fragend = int(row[4]) 
	desc = row[8]
	acc = extract_acc(desc)
	start = max(1, fragstart-flank) 
	end = min(contig_lengths_dict[contig], fragend + flank)
	cmd = ["esl-sfetch", "-o", "query.fna", "-n", desc, "-c", f"{start}..{end}", genome, contig]
	subprocess.run(cmd, check=True)
	query_outfile.write(desc + '\t' + contig + '\t' + str(start) + '\t' + str(end) + '\n')
	mainloc = main_pos_dict.get(acc, [])
	if len(mainloc) == 0: 
		continue
	os.system('rm target.fna')
	for exon in mainloc: 
		mainstart = int(exon[0]) 
		mainend = int(exon[1]) 
		maincontig = exon[2] 
		pos = exon[3]
		start = max(1, mainstart-flank)
		end = min(contig_lengths_dict[maincontig], mainend + flank) 
		numexons += 1
		outname = acc + '_' + pos 
		print(outname)
		eslcmd = 'esl-sfetch -n ' + outname + ' -c ' + str(start) + '..' + str(end) + ' ' + genome + ' ' + maincontig + ' >> target.fna'
		os.system(eslcmd)
		target_outfile.write(outname + '\t' + maincontig + '\t' + str(start) + '\t' + str(end) + '\n')
		#cmd = ["esl-sfetch", "-o", "target.fna", "-n", acc, "-c", f"{start}..{end}", genome, maincontig]
		#subprocess.run(cmd, check=True)
	dbcmd = 'makeblastdb -in target.fna -dbtype nucl'
	os.system(dbcmd)
	blastcmd = 'blastn -evalue 0.01 -dust no -task dc-megablast -query query.fna -db target.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length" >> ' + sys.argv[1] + '_source_blastout'
	os.system(blastcmd)
	os.system('rm target_CTRL.fna')
	for j in range(0, numexons):  # for every exon, pull a random other exon and do the same control
		randmain = random.choice(main_pos_list)
		randstart, randstop, ctrlcontig, ctrlacc = randmain
		ctrlstart = max(1, randstart-flank) 
		ctrlstop = min(contig_lengths_dict[ctrlcontig], randstop + flank)
		eslcmd = 'esl-sfetch -c ' + str(ctrlstart) + '..' + str(ctrlstop) + ' ' + genome + ' ' + ctrlcontig + ' >> target_CTRL.fna'
		# default name in esl is automatically 'contig'/'start'-'end'; this is sufficient
		os.system(eslcmd)
		#ctrl_outfile.write(ctrlacc + '\t' + ctrlcontig + '\t' + str(ctrlstart) + '\t' + str(ctrlstop) + '\n')
	dbcmd = 'makeblastdb -in target_CTRL.fna -dbtype nucl'
	os.system(dbcmd)
	blastcmd = 'blastn -evalue 0.01 -dust no -task dc-megablast -query query.fna -db target_CTRL.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length" >> ' + sys.argv[1] + '_CTRL_source_blastout'
	os.system(blastcmd) 
query_outfile.close()
target_outfile.close() 
#ctrl_outfile.close() 
