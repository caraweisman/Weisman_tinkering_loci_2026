import numpy as np 
import sys 
import re

consfile = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t', skip_header=True)

# dict that lists conservation status of each frag
# but formatting is fucked up relative to the frag listings here! so have to change first
# cons file eg: uncharacterized-protein-Dmel-CG5355-[Drosophila-melanogaster]_JAWNOL010000058.1:2375548-2376126
# orf file eg: uncharacterized-protein-Dmel-CG5355-[Drosophila-melanogaster]_2375548_2376126
frag_cons_dict = {} 
for row in consfile: 
	fragid = row[0]
	namestr = re.split('_', fragid)[0] # up to closing ]
	coords = re.split(':', fragid)[1] # start-end
	start = re.split('-', coords)[0]
	stop = re.split('-', coords)[1]
	newname1 = namestr + '_' + start + '_' + stop
	newname2 = namestr + '_' + stop + '_' + start
	cons_values = row[1:]
	frag_cons_dict[newname1] = cons_values
	frag_cons_dict[newname2] = cons_values # append both directions for start/stop since somehwo they got switched in the conservation ipipeline

outfile = open(sys.argv[1] + '_DAFF_SPECIFIC', 'w')
lenoutfile = open(sys.argv[1] + '_DAFF_SPECIFIC_ORFLENGTHS', 'w')

with open(sys.argv[1]) as f:
	codingorfs = [line.rstrip('\n').split('\t') for line in f]
with open(sys.argv[3]) as f:
	prot_length_lines = [line.rstrip('\n') for line in f]

totalorfs = 0
neworfs = 0
for k in range(0, len(codingorfs)): 
	row = codingorfs[k]
	print(row)
	# rows have different numbers of columns; first scan to see where all the fragments with positional information stop
	fragliststopidx = None
	for i in range(0, len(row)): 
		cell = row[i]
		if 'ORF' in cell: 
			fragliststopidx = i
			prefix=row[0:i] # prefix is information about the class of transcripts: what features they have and the total number of reads. 
			break 
	for i in range(0, len(row)):  
		cell = row[i]
		if 'ORF' in cell: 
			totalorfs += 1
			print(cell)
			cons = True
			complist = re.split(',',row[i+1])
			print(complist)
			for comp in complist: 
				print(comp)
				component_conservation = False # all fragments with the same name are assumed to share an origin; if one is conserved, all are conserved
				melstripped = re.split('\\[', comp)[0].replace('_','-') # weird formatting with dash before the bracket
				print(melstripped)
				for j in range(0, fragliststopidx): 
					if melstripped in row[j]: # frag entry with coordinates
						print(row[j])
						cons_array = frag_cons_dict.get(row[j], [])
						print('cons_array:', cons_array)
						if len(cons_array) == 0: 
							component_conservation = True # not in gff; pruned due to high e-value; assume true
						for species in cons_array: 
							if species == 'True': ## only looking for affinis-specific things!  
								component_conservation = True
				if component_conservation == False: 
					cons = False
			if cons == False: 
				neworfs += 1
				orf_entries = row[i:i+3] # i is the ORF call; i+1 is the components ; i +2 is the number of supporting reads
				relevant_entries = prefix + orf_entries # write first information about the transcript class, then information about this specific orf; per orf
				outfile.write("\t".join(relevant_entries))
				outfile.write('\n') # writes one line per ORF, not per read, like before; but lists all of the frags encoded at the nt level in that read, as before
				fields = prot_length_lines[totalorfs - 1].split('\t')
				if fields[0] != row[i+1]: 
					print('MISALIGN', totalorfs, fields[0], '!=', row[i+1])
				lenoutfile.write(prot_length_lines[totalorfs - 1] + '\n')


outfile.close()
lenoutfile.close()
print('TOTAL ORFS', totalorfs)
print('AFFINIS-SPECIFIC ORFS', neworfs)
