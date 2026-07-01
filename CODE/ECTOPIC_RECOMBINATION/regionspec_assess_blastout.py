import numpy as np 
import sys
import re 

gff =  np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)
blastout = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t', comments=None)
flank = int(sys.argv[3])
frag_conversion = np.genfromtxt(sys.argv[4], dtype=str, delimiter='\t', comments=None)
source_conversion = np.genfromtxt(sys.argv[5], dtype=str, delimiter='\t', comments=None)
control_blastout = np.genfromtxt(sys.argv[6], dtype=str, delimiter='\t', comments=None)
#control_source_conversion = np.genfromtxt(sys.argv[7], dtype=str, delimiter='\t', comments=None)
# not necessary! because I encoded the absolute start/stop coordinates in the subject ID for the control blast.

minalilen = 100

# in this script, enforce that the regions of similarity are around the same region of the protein
def get_frag_coords(desc): 
	coords = re.split('/', desc)[0][4:]
	start = int(re.split('-', coords)[0])
	stop = int(re.split('-', coords)[1])
	return start, stop

def get_source_coords(desc): # for matching exons between fragment/source
	coords = re.split('_', desc)[2]
	start = int(re.split('-', coords)[0])
	stop = int(re.split('-', coords)[1])
	return start, stop

#  dict to convert between fragment sequence and absolute coordinates/contig in genome
frag_conversion_dict = {} 
for row in frag_conversion: 
	fragid = row[0]
	contig = row[1]
	start = int(row[2]) 
	end = int(row[3])
	frag_conversion_dict[fragid] = [contig, start, end]

source_conversion_dict = {} 
for row in source_conversion: 
	sourceid = row[0]
	contig = row[1]
	start = int(row[2]) 
	end = int(row[3])
	source_conversion_dict[sourceid] = [contig, start, end]

blastdict = {} # format: key = fragment desc, value = list of blast hits 
for row in blastout: 
	fragid = row[0]
	if fragid not in blastdict: 
		blastdict[fragid] = []
	blastdict[fragid].append(row)

control_blastdict = {} #format: key = fragment desc, value = list of blast hits 
for row in control_blastout: 
	fragid = row[0]
	if fragid not in control_blastdict: 
		control_blastdict[fragid] = []
	control_blastdict[fragid].append(row)

outfile = open('recombination_hits', 'w')
outfile.write('#FRAG_ID \t FRAG_CONTIG \t FRAG_RECOMB_START \t FRAG_RECOMB_END \t SOURCE_ID \t SOURCE_CONTIG \t SOURCE_RECOMB_START \t SOURCE_RECOMB_END \n')

outfile2 = open('CONTROL_recombination_hits', 'w')
outfile2.write('#FRAG_ID \t FRAG_CONTIG \t FRAG_RECOMB_START \t FRAG_RECOMB_END \t SOURCE_CONTIG \t SOURCE_RECOMB_START \t SOURCE_RECOMB_END \n')

for row in gff: 
	fragid = row[8]
	frag_prot_start, frag_prot_stop = get_frag_coords(fragid) # protein range of fragment for finding the right source exon; pulls from query name, in which is embedded the coordinates of the fragment
	frag_contig, abs_frag_start, abs_frag_end  = frag_conversion_dict[fragid] # absolute start/end coordinates of the region around the fragment used in the blast search (separate file), because blast hit is relative to query
	blasthits = blastdict.get(fragid, [])
	for hit in blasthits: 
		sourceid = hit[4]
		source_prot_start, source_prot_stop = get_source_coords(sourceid) # protein range of source 
		source_contig, abs_source_start, abs_source_end  = source_conversion_dict[sourceid] # absolute start/end coordinates of source region, again bc blast coords are relative to nt query
		qlen = int(hit[1])
		# only look for matches that aren't in the fragment itself, ie are the first and last FLANK of the query, where FLANK is the size of the flanking region, chosen in bash script and passed here
		qprefix = [0, flank]
		qsuffix = [qlen-flank, qlen]
		qstart = int(hit[2])
		qend = int(hit[3])
		slen = int(hit[5]) 
		alilen = int(hit[10])
		sstart = min(int(hit[6]), int(hit[7]))
		send = max(int(hit[6]), int(hit[7]))
		sprefix = [0, flank]
		ssuffix = [slen-flank, slen]
		overlap = max(0, min(frag_prot_stop, source_prot_stop) - max(frag_prot_start, source_prot_start)) # make the fragment overlap the source exon by at least 1 aa (permissive but low null p)
		querypos = False
		subjpos = False
		if qstart < qprefix[1]: 
			effstop = min(qprefix[1], qend) 
			efflen = effstop - qstart 
			if efflen > minalilen: 
				querypos = True
		if qend > qsuffix[0]: 
			effstart = max(qsuffix[0], qstart) 
			efflen = qend - effstart
			if efflen > minalilen: 
				querypos = True
		if sstart < sprefix[1]:
			effstop = min(sprefix[1], send) 
			efflen = effstop - sstart
			if efflen > minalilen: 
				subjpos = True
		if send > ssuffix[0]: 
			effstart = max(ssuffix[0], sstart)
			efflen = send - effstart
			if efflen > minalilen: 
				subjpos = True
		if querypos == True and subjpos == True and overlap > 0: # fragment must come from this part of the protein
		#if qend < qprefix[1] and sstart < send and send < sprefix[1]: # both strands are plus relative to the features by design
			abs_fraghit_start = abs_frag_start + qstart # genomic start coordinate plus the relative start coordinates of the blast hit = genomic start coord of the blast hit
			abs_fraghit_end = abs_frag_start + qend # genomic end coordinate plus the relative start coordinates of the blast hit = genomic end coord of the blast hit
			abs_sourcehit_start = abs_source_start + sstart # same as two lines above for source fragment
			abs_sourcehit_end = abs_source_start + send
			outfile.write(fragid + '\t' + frag_contig + '\t' + str(abs_fraghit_start) + '\t' + str(abs_fraghit_end) + '\t' + sourceid + '\t' + source_contig + '\t' + str(abs_sourcehit_start) + '\t' + str(abs_sourcehit_end) + '\n')
		#elif qstart > qsuffix[0] and sstart < send and sstart > ssuffix[0]: 
		#	abs_fraghit_start = abs_frag_start + qstart
		#	abs_fraghit_end = abs_frag_start + qend
		#	abs_sourcehit_start = abs_source_start + sstart
		#	abs_sourcehit_end = abs_source_start + send
		#	outfile.write(fragid + '\t' + frag_contig + '\t' + str(abs_fraghit_start) + '\t' + str(abs_fraghit_end) + '\t' + source_contig + '\t' + str(abs_sourcehit_start) + '\t' + str(abs_sourcehit_end) + '\n')
	control_blasthits = control_blastdict.get(fragid, [])
	for hit in control_blasthits: 
		sourceid = hit[4]
		ctrl_source_contig = re.split('/', sourceid)[0]
		ctrl_abs_source_coords = re.split('/', sourceid)[1]
		ctrl_abs_source_start = int(re.split('-', ctrl_abs_source_coords)[0])
		ctrl_abs_source_stop = int(re.split('-', ctrl_abs_source_coords)[1]) # no conversion files necessary here because it's encoded into subject name in blast
		qlen = int(hit[1])
		qstart = int(hit[2])
		qend = int(hit[3])
		qsuffix = [qlen-flank, qlen]
		qprefix = [0, flank]
		slen = int(hit[5]) 
		alilen = int(hit[10])
		sstart = min(int(hit[6]), int(hit[7]))
		send = max(int(hit[6]), int(hit[7]))
		sprefix = [0, flank]
		ssuffix = [slen-flank, slen]
		querypos = False
		subjpos = False
		if qstart < qprefix[1]: 
			effstop = min(qprefix[1], qend) 
			efflen = effstop - qstart 
			if efflen > minalilen: 
				querypos = True
		if qend > qsuffix[0]: 
			effstart = max(qsuffix[0], qstart) 
			efflen = qend - effstart
			if efflen > minalilen: 
				querypos = True
		if sstart < sprefix[1]:
			effstop = min(sprefix[1], send) 
			efflen = effstop - sstart
			if efflen > minalilen: 
				subjpos = True
		if send > ssuffix[0]: 
			effstart = max(ssuffix[0], sstart)
			efflen = send - effstart
			if efflen > minalilen: 
				subjpos = True
		if querypos == True and subjpos == True:
			abs_fraghit_start = abs_frag_start + qstart
			abs_fraghit_end = abs_frag_start + qend  # these two are not really necessary to repeat, but makes script analysis easier / basic technical control
			abs_sourcehit_start = ctrl_abs_source_start + sstart
			abs_sourcehit_end = ctrl_abs_source_start + send
		#if qend < qprefix[1] and sstart < send and send < sprefix[1]: # both strands are plus relative to the features by design
			outfile2.write(fragid + '\t' + frag_contig + '\t' + str(abs_fraghit_start) + '\t' + str(abs_fraghit_end) + '\t' + sourceid + '\t' + ctrl_source_contig + '\t' + str(abs_sourcehit_start) + '\t' + str(abs_sourcehit_end) + '\n')
		#elif qstart > qsuffix[0] and sstart < send and sstart > ssuffix[0]: 
		#	abs_fraghit_start = abs_frag_start + qstart
		#	abs_fraghit_end = abs_frag_start + qend
		#	abs_sourcehit_start = abs_source_start + sstart
		#	abs_sourcehit_end = abs_source_start + send
		#	outfile.write(fragid + '\t' + frag_contig + '\t' + str(abs_fraghit_start) + '\t' + str(abs_fraghit_end) + '\t' + source_contig + '\t' + str(abs_sourcehit_start) + '\t' + str(abs_sourcehit_end) + '\n')
outfile.close() 
outfile2.close()
