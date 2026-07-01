import numpy as np 
import sys
import re

## usage: annot_to_gff.py (main hits file) (extra hits file) (outfile name)

main_hits_raw = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t')
extra_hits_raw = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t') 

id_col = np.full((main_hits_raw.shape[0], 1), 'main')
main_hits = np.hstack((id_col, main_hits_raw))

id_col = np.full((extra_hits_raw.shape[0], 1), 'extra')
extra_hits = np.hstack((id_col, extra_hits_raw))

full_annots_raw = np.vstack((main_hits, extra_hits))

# Convert column 8 (index 7) to numeric for sorting
start_col = full_annots_raw[:, 7].astype(float)

# Get sorted indices
#sorted_indices = np.argsort(start_col)
sorted_indices = np.lexsort((full_annots_raw[:, 7].astype(float), full_annots_raw[:, 5]))

# Apply sort to the full array
full_annots = full_annots_raw[sorted_indices]

contig_col = 5
start_col = 7
end_col = 8
prot_start_col = 3
prot_end_col = 4
prot_len_col = 2
prot_id_col = 1
eval_col = 9
strand_col = 12

outfile = open(sys.argv[3], 'w')

# shades of purple
#bacterial_colors = ['#BF40BF', '#5D3FD3', '#CF9FFF', '#DA70D6', '#E0B0FF']
# shades of gray
#TE_colors = ['#36454F', '#A9A9A9','#808080', '#D3D3D3']
#ROYGBPink, dark
#main_colors = ['#D2042D', '#FF5F1F', '#FFEA00','#7CFC00','#1F51FF','#FF00FF']
#ROYGBPink, pastel
#extra_colors = ['#FAA0A0', '#FAC898','#FFFAA0','#ECFFDC','#B6D0E2','#F8C8DC'] 

# shades of brown
bacterial_colors = ['#CD7F32','#C19A6B','#966919','#A0522D']
# shades of gray
TE_colors = ['#36454F', '#A9A9A9','#808080', '#D3D3D3']
# warm colors
extra_colors = ['#FF2400','#E97451', '#FFAC1C','#FDDA0D']
# cool colors
main_colors = ['#50C878', '#0BDA51', '#40E0D0', '#0096FF', '#5D3FD3', '#CF9FFF']



bactcolorpos = 0
TEcolorpos = 0
maincolorpos = 0
extracolorpos = 0

bactcoloridx = 0
TEcoloridx = 0
maincoloridx = 0
extracoloridx = 0

for i in range (0, len(full_annots)):
	entry = full_annots[i] 
	newrow = []
	newrow.append(entry[contig_col])
	newrow.append('tblastn')
	newrow.append('match')
	start = min(int(entry[start_col]), int(entry[end_col]))
	end = max(int(entry[start_col]), int(entry[end_col]))
	newrow.append(str(start))
	newrow.append(str(end))
	newrow.append('.')
	frame = int(re.split('/',entry[strand_col])[1])
	if frame > 0: 
		strand = '+'
	elif frame <0: 
		strand = '-'
	newrow.append(strand)
	newrow.append('.')
	namestr = '(' + entry[prot_start_col] + '-' + entry[prot_end_col] + '/' + entry[prot_len_col] + ')' + ' ' + entry[prot_id_col] + '_' + entry[start_col] + '_' + entry[end_col]
	descstr = 'evalue_' + entry[eval_col]
	previous_id = full_annots[i-1][prot_id_col] 
	isoformfree_previousid = re.split('isoform',re.split(' ', previous_id)[1])[0]
	curr_id = entry[prot_id_col]
	isoformfree_currid = re.split('isoform',re.split(' ', curr_id)[1])[0]
	prevframe = int(re.split('/',full_annots[i-1][strand_col])[1])
	if prevframe > 0:
		prevstrand = '+'
	else: 
		prevstrand = '-' 
	isoforms = False
	if isoformfree_previousid == isoformfree_currid and prevstrand == strand: 
		isoforms = True
	if entry[0] == 'main':
		if isoforms == True: 
			color = main_colors[maincoloridx]
		else: 
			maincolorpos += 1
			maincoloridx = maincolorpos % len(main_colors)
			color = main_colors[maincoloridx] 
	elif entry[0] == 'extra': 
		if 'DF' in entry[prot_id_col]: 
			if isoforms == True:
				color = TE_colors[TEcoloridx]
			else: 
				TEcolorpos += 1
				TEcoloridx = TEcolorpos % len(TE_colors) 
				color = TE_colors[TEcoloridx]
		elif 'WP' in entry[prot_id_col]: 
			if isoforms == True: 
				color = bacterial_colors[bactcoloridx]
			else: 
				bactcolorpos += 1
				bactcoloridx = bactcolorpos % len(bacterial_colors) 
				color = bacterial_colors[bactcoloridx]
		else: 
			if isoforms == True: 
				color = extra_colors[extracoloridx]
			else: 
				extracolorpos += 1
				extracoloridx = extracolorpos % len(extra_colors)
				color = extra_colors[extracoloridx]
	fullfield = 'ID=' + namestr + ';Name=' + namestr + ';description=' + descstr + ';color=' + color
	encoded_fullfield = fullfield.replace(' ', '-')
	newrow.append(encoded_fullfield) 
	outfile.write('\t'.join(newrow) + '\n')
outfile.close()
