import sys
import os
import re

# usage = find_reads_by_cons.py (samtools-derived gff) (fasta) 
# takes 31 cores

fasta = sys.argv[2]

def extract_name(desc):
    ID=re.split(';', desc)[0]
    intacc = re.split('[-_]', ID)[4:-2]
    if 'isoform' in intacc: 
        isoidx = intacc.index('isoform')
        intacc.pop(isoidx+1)
        intacc.pop(isoidx)
        intacc[isoidx-1] = intacc[isoidx-1][:-1]
    name = '-'.join(intacc)
    if 'None' in name: 
        print(desc)
    return name

def extract_acc(desc): 
	acc = re.split('-', desc)[2]
	return acc


prot_reads_dict = {}
with open(sys.argv[1]) as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		read = fields[3]
		desc = fields[20]
		acc = extract_acc(desc) 
		if acc not in prot_reads_dict: 
			prot_reads_dict[acc] = set()
		prot_reads_dict[acc].add(read)

prot_reads_dict = {acc: list(reads) for acc, reads in prot_reads_dict.items()}

#cmd = 'mkdir REP_READS'
#os.system(cmd)

dbcmd = 'esl-sfetch --index ' + fasta
os.system(dbcmd)

# took forever to run and timed out so restarting based on where it left off
with open('done_accs.txt') as f:
    done_list = set(line.strip() for line in f if line.strip())

for prot in prot_reads_dict: 
	if prot in done_list: 
		continue
	clearcmd = 'rm -rf tmp'
	os.system(clearcmd) # clear tmp files from previous run
	reads = prot_reads_dict[prot]
	outfile = open('readlist', 'w')
	for read in reads: 
		outfile.write(read + '\n') 
	outfile.close() 
	eslcmd = 'esl-sfetch -f ' + fasta + ' readlist > readlist.fna'
	os.system(eslcmd)
	# extract all seqs overlapping a protein
	clustcmd = 'mmseqs easy-cluster readlist.fna ' + prot + ' tmp --min-seq-id 0.95 -c 0.3 --cov-mode 1 --cluster-mode 2 --threads 31'
	os.system(clustcmd)
	# cluster
	mvcmd = 'mv ' + prot + '*' + ' REP_READS_2/'
	os.system(mvcmd)
	# move to subdir for cleanliness
	
