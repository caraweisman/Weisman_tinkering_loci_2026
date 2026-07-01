import numpy as np 
import sys
import re 

unsorted_gff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)

# sort first by read, then by contig, then by start position of read on feature 
# i think it's already sorted this way, but just to be sure
gff = unsorted_gff[np.lexsort((unsorted_gff[:, 1].astype(int), unsorted_gff[:, 0], unsorted_gff[:, 3]))]

polyA_ts_strands = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t') # backup ts strands (transcript orientation) derived from polyA tails for cases where minimap can't deduce from splice

# put strands derived from polyA polarization into dict for easy lookup
polyA_ts_strands_dict = {}
for row in polyA_ts_strands: 
	read = row[0]
	strand = row[1]
	polyA_ts_strands_dict[read] = strand

main_colors = ['#50C878', '#0BDA51', '#40E0D0', '#0096FF', '#5D3FD3', '#CF9FFF']
frag_colors = ['#FF2400','#E97451', '#FFAC1C','#FDDA0D']

def extract_name(desc):
	ID=re.split(';', desc)[0]
	tokens = re.split('[-_]', ID)
	intacc = tokens[4:-2]
	if 'isoform' in intacc: 
		isoidx = intacc.index('isoform')
		intacc.pop(isoidx+1)
		intacc.pop(isoidx)
		intacc[isoidx-1] = intacc[isoidx-1][:-1]
	name = '-'.join(intacc) + '_' + '_'.join(tokens[-2:])
	if 'None' in name: 
		print(desc)
	return name

def extract_acc(desc): 
	acc = re.split('-', desc)[2]
	return acc


clusters = {}
with open('clusters.txt') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		focal = fields[0]
		targets = fields[1:]
		clusters[focal] = targets


read_features = {} 
# for each read: [[name, start, end, strand, type (cons/frag)]...]


# for each read, make a list of the features that it overlaps, in order of their appearance on the read
# for now, use the order in the gff thus assuming all reads are plus, but store strand to later reorient from 5' to 3' on each read
# deduplicate overlapping features for conserved genes to avoid multiple isoforms being appended 
for row in gff: 
	read = row[3]
	if read not in read_features: 
		read_features[read] = [] # initialize lists for main and fragment components
	featstart = int(row[1]) 
	featend = int(row[2]) # not really the feature but the read's overlap with it
	featlen = featend-featstart
	#strand = row[5] # + or -; strand of READ, not feature 
	strand = row[21] # "bio": what i want: biological oreintation of read on reference. derived from operation on ts/fs: if same, +; if different, -
	if strand == '.':  # if this isn't assigned, it's because minimap couldn't determine read orientation because no sufficient splie sites; use polarization from polyA tail in sys.argv[2]
		fs = row[22] # fs is always available from samtools; i reprinted it here
		ts = polyA_ts_strands_dict[read] # get ts from polyA calls
		if fs == ts: 
			strand = '+'
		else: 
			strand = '-'
	desc = row[20]
	name = extract_name(desc)
	acc = extract_acc(desc) 
	if any(s in desc for s in main_colors) == True:
		info = [name, featstart, featend, strand, 'CONS', acc]
		overlap = False
		# deduplicate with other elements to avoid a gazillion isoforms of the same thing being appended and so on - should mostly affect conserved genes for this reason
		for element in read_features[read]: # for other appended main elements that the read overlaps
			element_start = element[1]
			element_end = element[2]
			element_type = element[4]
			element_len = element_end - element_start
			overlaplen = max(0, min(element_end, featend) - max(element_start, featstart))
			if float(overlaplen)/featlen > 0.5 or float(overlaplen)/element_len > 0.5: # doesn't matter here whether overlap is between cons/cons or cons/frag
				overlap = True
		if overlap == False:
			read_features[read].append(info)
	elif any(s in desc for s in frag_colors) == True:
		info = [name, featstart, featend, strand, 'FRAG', acc]
		read_features[read].append(info)



feature_groups = {} # dictionary of 5' to 3' on read ordered elements
# keys are an ordered list of the form [[element1, type(cons/frag)], [element2, type], ...]
# values are read ids
# based on ordered gff, elements should appear in their order on the contig; have to correct for read strand
for read in read_features:
	ordered_features = []
	featureinfos = read_features[read]
	strand = featureinfos[0][3] # property of the read, so will be the same for all of a read's features; pick the first
	print(read, repr(strand), len(featureinfos), file=sys.stderr)   #### MOVED + EXPANDED
	print(read, repr(strand), len(featureinfos), 'OFI' if 'ordered_featureinfos' in dir() else 'UNSET', file=sys.stderr)
	if strand == '+': 
		ordered_featureinfos = featureinfos
	elif strand == '-':
		ordered_featureinfos = featureinfos[::-1]
	else:
		sys.stderr.write(f"skipping read {read}: strand still neither -/+ after polyA fill \n")
		continue
	for elem in ordered_featureinfos: 
		name = elem[0]
		elemtype = elem[4]
		acc = elem[5]
		ordered_features.append(name)
		ordered_features.append(elemtype)
		ordered_features.append(acc) # for purposes of clustering
	featureskey = tuple(ordered_features)
	if featureskey not in feature_groups: 
		feature_groups[featureskey] = [read]
	else: 
		feature_groups[featureskey].append(read) 

outfile = open(sys.argv[1] + '_analysis', 'w') 
outfile2 = open(sys.argv[1] + '_analysis_withreads', 'w') 

outfile.write('#READ_TYPE \t CONS_NUM \t FRAG_NUM \t NUM_READS \t IDS \n')
outfile2.write('#READ_TYPE \t CONS_NUM \t FRAG_NUM \t IDS \t READS \n')

for feature in feature_groups: 
	consnum = 0
	fragnum = 0
	conselements = []
	consaccs = []
	fragelements = []
	fragaccs = []
	for i in range(0, len(feature)-2): # groups of 3; exclude last two
		if i % 3 == 0: # every third  elements are the names; odd are the types
			element = feature[i]
			elemtype = feature[i+1]
			acc = feature[i+2]
			if elemtype == 'CONS': 
				if element not in conselements:
					consnum += 1
					conselements.append(element)
					consaccs.append(acc)
			if elemtype == 'FRAG': 
				if element not in fragelements:
					fragnum += 1
					fragelements.append(element)
					fragaccs.append(acc)
	if consnum == 0 and fragnum == 1: 
		readtype = 'SINGLETON'
	elif consnum == 0 and fragnum > 1:
		purgedfraglist = [] # check to see if the fragments are in the same cluster; if so, purge from list
		purgedacclist = []
		for j in range(0, len(fragelements)):
			frag = fragelements[j]
			acc = fragaccs[j] 
			inclust = False
			for k in range(0, len(purgedfraglist)):
				prevacc = purgedacclist[k]
				clust = clusters[prevacc]
				if acc in clust:  
					inclust = True
			if inclust == False: 
				purgedfraglist.append(frag)
				purgedacclist.append(acc)
		fragnum = len(purgedfraglist) # recalculate fragnum
		if fragnum > 1: 
			readtype = 'CHIMERA' 
		elif fragnum == 1: 
			readtype = 'SINGLETON'
	elif consnum > 0 and fragnum > 0:
		if len(set(fragelements))<=1 and len(set(conselements))<=1 and conselements[0] == fragelements[0]: # remove cases where it's all the same gene because of annotation fuzz
			continue # don't even bother writing these
		purgedfraglist = []
		for j in range(0, len(fragelements)): 
			frag = fragelements[j]
			acc = fragaccs[j] 
			inclust = False
			for k in range(0, len(conselements)): 
				consacc = consaccs[k]
				clust = clusters[consacc]
				if acc in clust: 
					inclust = True
			if inclust == False: 
				purgedfraglist.append(frag) 
		fragnum = len(purgedfraglist) 
		if fragnum > 0: 
			readtype = 'NEW_EXON'
		elif fragnum == 0:
			continue
	else: 
		readtype = 'OTHER'
	outfile.write(readtype + '\t' + str(consnum) + '\t' + str(fragnum))
	outfile2.write(readtype + '\t' + str(consnum) + '\t' + str(fragnum))
	reads = feature_groups[feature]
	outfile.write('\t' + str(len(reads)))
	for z in range(0, len(feature)):
		element = feature[z] 
		if z % 3 == 0 or z % 3 == 1 or z % 3 == 2: # previously just had first two ifs to prevent accession being printed but now it's useful to include and too lazy to dedent
			outfile.write('\t' + element)
			outfile2.write('\t' + element)
	readlist = ','.join(reads)
	outfile2.write('\t' + readlist + '\n')
	outfile.write('\n')

outfile.close()
outfile2.close()
