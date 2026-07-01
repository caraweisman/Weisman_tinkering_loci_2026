import numpy as np 
import sys
import os
import re 
import subprocess
import random

fragff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)
flank = int(sys.argv[2])
contig_lengths = np.genfromtxt(sys.argv[3], dtype=str)
genome = sys.argv[4]

contig_lengths_dict = {}
for contig in contig_lengths: 
    contigname = contig[0]
    contiglen = int(contig[1])
    contig_lengths_dict[contigname] = contiglen

tx_paths = ['/n/eddy_lab/users/cweisman/TINKERING_LOCI/G0_2V3/AFFINIS/ANALYSIS/CONS_TRANSCRIPTION/REP_READS', '/n/eddy_lab/users/cweisman/TINKERING_LOCI/G0_2V3/AFFINIS/ANALYSIS/CONS_TRANSCRIPTION/REP_READS_2']

def extract_acc(desc): 
	acc = re.split('-', desc)[2]
	return acc

def extract_pos(desc):
	match = re.search(r'\((\d+-\d+)/', desc)
	result = match.group(1)  # '1-652'
	return result 


os.system('rm ' + sys.argv[1] + '_TX_blastout') # clear old blast outputs
#os.system('rm ' + sys.argv[1] + '_CTRL_source_blastout') # clear old blast outputs
#query_outfile = open(sys.argv[1] + '_fragment_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later  
#target_outfile = open(sys.argv[1] + '_source_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later  
#ctrl_outfile =  open(sys.argv[1] + '_CTRL_source_blastout_coords', 'w') # need a file to keep track of the coordinates of the searched regions to look for TEs/triggers later 
#not actually necessary because coordinates are included in target id; see next script

# take every fragemnt, name it according to its description in the gff, add 2500 nt of flanking sequence (if contig length allows), and blast it against transcript from parental gene (isoform)
tx_not_found_list = []
for i in range(0, len(fragff)):
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
	#query_outfile.write(desc + '\t' + contig + '\t' + str(start) + '\t' + str(end) + '\n')
	putdbname = acc + '_rep_seq.fasta' # name of preassembled dbs of repreesentative transcripts per isoform from gff and rnaseqdata
	found = False
	for path in tx_paths: 
		filepath = path + '/' + putdbname
		if os.path.exists(filepath) == True: 
			cmd = 'cp ' + filepath + ' target.fna'
			os.system(cmd)
			found = True
	if found == False: 
		print('ACC TX NOT FOUND:', acc)
		tx_not_found_list.append(desc)
		continue
	dbcmd = 'makeblastdb -in target.fna -dbtype nucl'
	os.system(dbcmd)
	blastcmd = 'blastn -evalue 0.01 -dust no -task dc-megablast -query query.fna -db target.fna -outfmt="6 qseqid qlen qstart qend sseqid slen sstart send evalue pident length" >> ' + sys.argv[1] + '_TX_blastout'
	os.system(blastcmd)
#query_outfile.close()
#target_outfile.close() 
#ctrl_outfile.close() 

outfile = open('no_tx_found', 'w')
for desc in tx_not_found_list: 
	outfile.write(desc + '\n')
outfile.close() 
