import numpy as np
import re 

# this version doesn't exclude coding regions from the transcript other than the matching exon (But does still only look for hits matching the flanking regions of framgnet)

blastout=np.genfromtxt('TL_frags.gff_TX_blastout', dtype=str, delimiter='\t', comments=None)
TLs = np.genfromtxt('TL_frags.gff', dtype=str, delimiter='\t', comments=None)

# make a dictionary of all noncoding regions from each representative transcript
# this will answer the question of whether the DNA-level hits from the recombination search are equally well explained by noncoding exons of the transcript
# coding exon that is the source of the fragment has already been removed in the dna search; also rmove here



def get_source_coords(desc): # for matching exons between fragment/source. one isoform per read and one length per isoform so just need start/end
	coords = re.split('/', desc)[0][4:]
	start = int(re.split('-', coords)[0])
	stop = int(re.split('-', coords)[1])
	return start, stop


read_codingregion_dict = {} # key = readid, value = [[readposstart, readposstop, exonprotstart, exonprotstop]..]: array of bounds of coding regions, to be selected and then avoided 
with open('feature_positions_on_reads_REPREADS.tsv') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		read = fields[0]
		featstart = int(fields[4])
		featend = int(fields[5])
		desc = fields[8] # feature description 
		protstart, protstop = get_source_coords(desc)
		if read not in read_codingregion_dict: 
			read_codingregion_dict[read] = []
		read_codingregion_dict[read].append([featstart, featend, protstart, protstop])
		# don't put read length here because it isn't in this fucking tsv; take it from blast output

minalilen = 100

TL_dict = {}
for TL in TLs: 
	ID=TL[8]
	TL_dict[ID] = []

def complement(ranges, universe):
	lo, hi = universe
	result = []
	cursor = lo
	for a, b in sorted(ranges):
		if a > cursor:
			result.append([cursor, a - 1])
		cursor = max(cursor, b + 1)
	if cursor <= hi:
		result.append([cursor, hi])
	return result

# file is named full_TX_hits; other version, which exlcudes coding regions, is named TX_hits
outfile = open('full_TX_hits', 'w')
for row in blastout: 
	fragid = row[0]
	fragprotstart, fragprotend = get_source_coords(fragid)
	hitstart = int(row[2])
	hitend = int(row[3]) 
	qlen = int(row[1]) 
	sstart = min(int(row[6]), int(row[7]))
	send = max(int(row[6]), int(row[7]))
	read = row[4]
	readlen = int(row[5])
	prefix = [0, 2500]
	suffix = [qlen-2500, qlen]
	hitlen = hitend-hitstart
	qrange = False
	srange = False
	if hitstart < prefix[1]: 
		realstop = min(hitend, prefix[1])
		efflen = realstop - hitstart
		if efflen > minalilen: 
			qrange = True
			#outfile.write('\t'.join(row) + '\n')
	if hitend > suffix[0]: 
		realstart = max(hitstart, suffix[0]) 
		efflen = hitend - realstart
		if efflen > minalilen: 
			qrange = True
			#outfile.write('\t'.join(row) + '\n')
	# only do read lookup if worthwile - expensive
	if qrange == True: # excludes hits just to fragment boundaries
		exonranges = read_codingregion_dict.get(read, [])
		allexcluded = []
		for i in range(0, len(exonranges)):
			readposstart, readposstop, exonprotstart, exonprotstop = exonranges[i]
			overlap = max(0, min(fragprotend, exonprotstop) - max(fragprotstart, exonprotstart)) # make the fragment overlap the source exon by at least 1 aa (permissive but low null p) 
			if overlap > 0: # correct exon
				excluded = [min(readposstart, readposstop), max(readposstart, readposstop)]
				allexcluded.append(excluded) 
		allowedranges = complement(allexcluded, [1, readlen]) # allowed range for overlap on read
		for arange in allowedranges: 
			overlap = max(0, min(arange[1], send) - max(arange[0], sstart) + 1)
			if overlap > minalilen:
				outfile.write(fragid + '\t' + str(hitstart) + '\t' + str(hitend) + '\t' + read + '\t' + str(sstart) + '\t' + str(send) + '\n')
				break # one is enough

outfile.close()

