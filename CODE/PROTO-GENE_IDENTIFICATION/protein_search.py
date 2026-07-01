import numpy as np
import sys
import os
import re

# updated 6/1 to print anotehr file with orf lengths: per configuration, list of max orf lengths per read for that configuration. to compare to Zhou 2008 which relies on 
# 2007 Flybase annotation which uses a cutoff of 50aa for annotated genes with est evidence and homology to tohe rproteins.

# usage: python (script.py) (readfile) (peptide file from esl-translate) (refdb) 

# readfile: file resulting from find_chimeric_reads.py with list of read type and fragments in it
# peptide file: esl-translate file of the reads
# refdb: melanogaster proteins

refdb = sys.argv[3]


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

def extract_name_blast(desc):
	# desc looks like: "NP_001303394.1 nicotinic acetylcholine receptor alpha4, isoform H [Drosophila melanogaster]"
	# strip the leading accession (first whitespace token)
	rest = desc.split(' ', 1)[1] if ' ' in desc else desc
	# remove ", isoform X" (X is one or more non-bracket, non-comma chars)
	rest = re.sub(r',\s*isoform\s+[^,\[]+', '', rest)
	# collapse all whitespace to single dashes
	name = re.sub(r'\s+', '-', rest.strip())
	return name


def extract_acc(desc): 
	acc = re.split(' ', desc)[0]
	return acc


# index for sfetch
eslcmd = 'esl-sfetch --index ' + refdb
os.system(eslcmd)
eslcmd = 'esl-sfetch --index ' + sys.argv[2]
os.system(eslcmd)

with open(sys.argv[2]) as f:
	peptidefile = [line.split() for line in f if line.startswith('>')]


# dictionary to convert from read id to a list of its esl-translate peptides
reads_to_eslpeps = {}
eslpeps_to_reads = {}
for row in peptidefile: 
	readid = re.split('=',row[1])[1]
	orfid = row[0][1:]
	if readid not in reads_to_eslpeps: 
		reads_to_eslpeps[readid] = []
	reads_to_eslpeps[readid].append(orfid) 
	if orfid not in eslpeps_to_reads: 
		eslpeps_to_reads[orfid] = []
	eslpeps_to_reads[orfid].append(readid)

finaloutfile = open(sys.argv[1] + '_CODINGORFS', 'w', buffering=1)

with open(sys.argv[1]) as f:
	readfile = [line.rstrip('\n').split('\t') for line in f if line.strip()][1:]

TESTreadoutfile = open('TEST_ORF_readlist_' + sys.argv[2], 'w')

chimericorfsearchinfo = open('CHIMERIC_ORF_POSITION_INFO_' + sys.argv[2], 'w', buffering=1)

lengthinfo = open('PROT_LENGTH_INFO_' + sys.argv[2], 'w', buffering=1)

for row in readfile: 
	classification = row[0]
	readnum = str(len(re.split(',',row[-1])))
	finaloutfile.write('\t'.join(row[0:-1])) # repeat current row so you knw where result is coming from except the read list which is unwieldy
	finaloutfile.write('\t' + readnum)
	TESTreadoutfile.write('\t'.join(row[0:-1])) # repeat current row so you knw where result is coming from except the read list which is unwieldy
	TESTreadoutfile.write('\t' + readnum)
	category = row[0]
	refprots = {} # dict: accession, type. list of all protein ids covered by the read and their type
	for i in range(0, len(row[0:-1])):
		field = row[i]
		if 'NP_' in field or 'XP_' in field or 'YP_' in field:
			if field not in refprots:
				refprots[field] = row[i-1] # previous entry is type, CONS/FRAG
	# make fasta of all proteins in the read for blast to see if they have an intact ORF in every read with these features
	outfile = open('REFACCS', 'w')
	for refprot in refprots: 
		outfile.write(refprot + '\n')
	outfile.close() 
	eslcmd = 'esl-sfetch -f ' + refdb + ' REFACCS > REFS.fa'
	os.system(eslcmd)
	print(eslcmd)
	readlist = re.split(',',row[-1])
	# make list of all intact ORFs in these reads from esl-translate
	orfs = []
	for read in readlist: 
		readorfs = reads_to_eslpeps.get(read, [])
		orfs = orfs + readorfs
	outfile = open('ORFIDS', 'w') 
	for orf in orfs: 
		outfile.write(orf + '\n') 
	outfile.close() 
	eslcmd = 'esl-sfetch -f ' + sys.argv[2] + ' ORFIDS > ORFS.fa'
	os.system(eslcmd) 
	print(eslcmd)
	dbcmd = 'makeblastdb -in REFS.fa -dbtype prot'
	os.system(dbcmd)
	# do blast search between proteins covered by reads and intact ORFs from reads
	blastcmd = 'blastp -query ORFS.fa -db REFS.fa -evalue 0.01 -outfmt="6 qseqid qlen qstart qend salltitles slen sstart send evalue pident length frames" > CHIMERA_BLASTOUT' 
	os.system(blastcmd) 
	if os.path.getsize('CHIMERA_BLASTOUT') == 0:
		finaloutfile.write('\n')
		continue
	blastout = np.genfromtxt('CHIMERA_BLASTOUT', dtype=str, delimiter='\t', ndmin=2)
	orfhits_dict = {} # dictionary for protein hits. format: key = orf, value = [protein hits]
	orf_blastdict = {} # dictionary for blast rows: format: key = orf, value = [[blast rows]]
	prots_name_to_acc = {} # dictionary to convert between protein hits' accession and name; both necessary in different places
	for blastrow in blastout: 
		orf = blastrow[0]
		hitname = extract_name_blast(blastrow[4]) # full info converted to isoform-free name
		hitacc = extract_acc(blastrow[4]) # full info converted to accession
		hitinfo = [hitname, hitacc]
		if hitname not in prots_name_to_acc: 
			prots_name_to_acc[hitname] = hitacc
		if orf not in orfhits_dict: 
			orfhits_dict[orf] = []
		if orf not in orf_blastdict: 
			orf_blastdict[orf] = []
		orf_blastdict[orf].append(blastrow)
		# keep track of which proteins are in which ORFs
		if hitname not in orfhits_dict[orf]:
			orfhits_dict[orf].append(hitname)
	seenconfigs = {} # convert to dictionary that lists different protein configurations. key = (hits) (tuple bc dict needs for key), value = [orfs with that config]
	seenconfig_lengths = {} # key = hits, value = [lengths of orfs]; similar to above for length info
	for orf in orfhits_dict:
		hits = tuple(sorted(orfhits_dict[orf])) # tuple of proteins found in ORF; sorted so they can be collapsed
		if hits not in seenconfigs: 
			seenconfigs[hits] = []
		if hits not in seenconfig_lengths: 
			seenconfig_lengths[hits] = []
		seenconfigs[hits].append(orf) 
		orflen = orf_blastdict[orf][0][1] # second column is qlen; choose first abtirarily since there may only be one
		seenconfig_lengths[hits].append(orflen)
		if len(hits) > 1 and classification == 'CHIMERA': # for chimeras, write all blast rows to output file for use in finding hcimeras in outgroup genomes
			# need and classification gate because classification takes clustering into account, which is not apparent from the list of fragments
			for blastrow in orf_blastdict[orf]: 
				chimericorfsearchinfo.write('\t'.join(blastrow) + '\n')
	# for every combination of proteins seen in this feature category
	for config in seenconfigs: 
		# for every orf with that combination, find its read and append it to totalreadlist, the list of reads with that configuratin (orf with constituent proteins)
		configorfs = seenconfigs[config]
		configstring = ','.join(config)
		totalreadlist = []
		for orf in configorfs: 
			readlist = eslpeps_to_reads.get(orf, [])
			for read in readlist: 
				if read not in totalreadlist: 
					totalreadlist.append(read)
		longest_per_read = {} # read -> (length, orfid); keep only each read's longest orf for this config
		for orf in configorfs: 
			orflen = int(orf_blastdict[orf][0][1])
			for read in eslpeps_to_reads.get(orf, []): 
				if read not in longest_per_read or orflen > longest_per_read[read][0]: 
					longest_per_read[read] = (orflen, orf)
		orflens = [str(longest_per_read[read][0]) for read in longest_per_read]
		prottypes = []
		for prot in config: 
			prottype = refprots[prots_name_to_acc[prot]] # dict takes accessio mto type
			prottypes.append(prottype) 
		if prottypes.count('FRAG') == 1 and prottypes.count('CONS') == 0:
			finaloutfile.write('\t' + 'ORF_SINGLETON' + '\t' + configstring + '\t' + str(len(totalreadlist))) 
			TESTreadoutfile.write('\t' + 'ORF_SINGLETON' + '\t' + configstring)
			for read in totalreadlist: 
				TESTreadoutfile.write('\t' + read)
			TESTreadoutfile.write('\n')
			lengthinfo.write(configstring + '\t' + '\t'.join(orflens) + '\n')
		elif prottypes.count('FRAG') > 1 and prottypes.count('CONS') == 0 and classification != 'SINGLETON': # classification has clustering niformation baked in that counting accesions doesn't; overrules
			finaloutfile.write('\t' + 'ORF_CHIMERA' + '\t' + configstring + '\t' + str(len(totalreadlist)))
			TESTreadoutfile.write('\t' + 'ORF_CHIMERA' + '\t' + configstring)
			for read in totalreadlist: 
				TESTreadoutfile.write('\t' + read)
			TESTreadoutfile.write('\n')
			lengthinfo.write(configstring + '\t' + '\t'.join(orflens) + '\n')
		elif prottypes.count('FRAG') > 0 and prottypes.count('CONS') > 0 and classification != 'SINGLETON': # classification has clustering niformation baked in that counting accesions doesn't; overrules
			finaloutfile.write('\t' + 'ORF_NEW_EXON' + '\t' + configstring + '\t' + str(len(totalreadlist))) 
			TESTreadoutfile.write('\t' + 'ORF_NEW_EXON' + '\t' + configstring)
			for read in totalreadlist: 
				TESTreadoutfile.write('\t' + read)
			TESTreadoutfile.write('\n')
			lengthinfo.write(configstring + '\t' + '\t'.join(orflens) + '\n')
	finaloutfile.write('\n') 
	finaloutfile.flush()
finaloutfile.close()
TESTreadoutfile.close()
chimericorfsearchinfo.close()
lengthinfo.close()

# final output in main outfile is: 
# status of feature list
# count of cons/frag features
# list of all features: frag status, accession, name
# total number of reads with those features 
# list of: ORF status / configuration / reads with that configuration; sum of reads per configuration should be equal or less than the total number of reads with those features on the tx
