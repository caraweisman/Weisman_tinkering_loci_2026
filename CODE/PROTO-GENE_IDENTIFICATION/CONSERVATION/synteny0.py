import numpy as np
import sys
import re
import math

# updated 5/16/26: changed to avoid errors when the fragment in affinis/focal is annotated as conserved in other gffs, but everything is otherwise the same.
# updated 5/17/26: changed from only looking for fragment in the case of two conserved flanks to adding one conserved flank cases and writing those regions to the target species esl output
# and: there are two cons flanks in focal but the target search still fails, fall through to the next test, where each flank is searched independently, on either side
# and: uses clustering. CLUSTER FILE MUST BE IN DIR; THIS IS HARD CODED 
# updated 5/20/26: if synteny not found because target synteny is not preserved but there are still two flanking genes, write both to output file, vs just the case where one flank exists; allows up to 2 writes to coords per frag

# finds flanking conserved genes for fragments in gff1
# usage: python (code.py) gff1 gff2 contigs

gff1 = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)
gff2 = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t', comments=None)

# sort gff first by contig and then by position!
gff1 = gff1[np.lexsort((gff1[:, 3].astype(np.int64), gff1[:, 0]))]
gff2 = gff2[np.lexsort((gff2[:, 3].astype(np.int64), gff2[:, 0]))]

main_colors = ['#50C878', '#0BDA51', '#40E0D0', '#0096FF', '#5D3FD3', '#CF9FFF']
frag_colors = ['#FF2400','#E97451', '#FFAC1C','#FDDA0D']

# initialize dictionaries for each gff: each contig gets an array, to be populated by arrays, 
# one for each conserved gene, listing the gff rows in which its fragments are found

contig_lengths_file1 = np.genfromtxt(sys.argv[3], dtype=str) # need contig lengths here to identify syntenic interval; 

clusterfile = 'clusters.txt' ### HARD CODED, MUST BE IN DIR

contig_lengths = {}
for contig in contig_lengths_file1: 
    contigname = contig[0]
    contiglen = int(contig[1]) 
    contig_lengths[contigname] = contiglen

gff1_contig_main_dict = {}
for contig in contig_lengths_file1: 
    contigname = contig[0]
    gff1_contig_main_dict[contigname] = [] 

gff2_contigs = np.unique(gff2[:, 0]) # don't need lengths, so just take all contigs in file
gff2_contig_main_dict = {}
for contig in gff2_contigs: 
    contigname = contig
    gff2_contig_main_dict[contigname] = [] 

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

# clusters are by accession
cluster_dict = {}
with open(clusterfile) as f: 
    for line in f:
        fields = line.rstrip('\n').split('\t')
        focprot = fields[0].strip()
        cluster_dict[focprot] = []
        for hit in fields[1:]: 
            cluster_dict[focprot].append(hit.strip())


# make a dictionary of the conserved genes in gff1
# key is gene name; value is an array of arrays, [contig, row], w multiple rows per multi-exon gene
for i in range(0, len(gff1)):
    row = gff1[i]
    contig = row[0]
    if any(s in row[8] for s in main_colors):
        protname = extract_name(row[8]) 
        protacc = extract_acc(row[8])
        gff1_contig_main_dict[contig].append([i, protname, protacc])

# make a dictionary of the conserved genes in gff2
# key is gene name; value is an array of arrays, [contig, row],  w multiple rows per multi-exon gene
for i in range(0, len(gff2)):
    row = gff2[i]
    contig = row[0]
    if any(s in row[8] for s in main_colors):
        protname = extract_name(row[8]) 
        protacc = extract_acc(row[8])
        gff2_contig_main_dict[contig].append([i, protname, protacc])


outfile1 = open('focal_fragment_synteny_eslcoords', 'w', buffering=1)
outfile2 = open('target_fragment_synteny_eslcoords', 'w', buffering=1)
outfile3 = open('synteny_broken_list', 'w', buffering=1)
outfile4 = open('0_syntenyfound_list', 'w', buffering=1)

# first,
# go through gff1 and for every fragment: 
# find the contig and coordinates of a "syntenic region", used for doing a blast search to find the fragment. 
# if the fragment is between two conserved exons, use their start (of upstream) and end (of downstream) - include them. 
# if the fragment is between a conserved exon and the start or end of a contig, 
# use the start/end of the conserved exon (to include it) and end/start of the contig as the coordinate.
# if the fragment is on its own contig with no conserved genes, use the start/end boundaries of that contig. (maybe terrible idea.)

totalrows = len(gff1)
for i in range(0, len(gff1)):
    if i % 100 == 0: 
        print(i, 'out of ', totalrows)
    synteny_found = False
    focal_flanks = False
    target_flanks_syntenic = False
    row = gff1[i]
    if any(s in row[8] for s in frag_colors) == True:
        startpos = None
        endpos = None
        fragrow = i
        print(row[8])
        fragid = extract_name(row[8])
        frag_acc = extract_acc(row[8])
        contig = row[0]
        fragstart = row[3]
        fragend = row[4]
        eslid = fragid + '_' + contig + ':' + fragstart + '-' + fragend
        contiglist = gff1_contig_main_dict[contig]
        leftcons = None
        rightcons = None
        for j in range(0, len(contiglist)): 
            consrow = contiglist[j][0]
            consrowstart = gff1[consrow][3]
            consrowend = gff1[consrow][4]
            protname = contiglist[j][1]
            protacc = contiglist[j][2]
            if j == 0 and consrow > fragrow: # if the first conserved gene on the contig is already farther than the fragment: 
                rightcons = protname
                rightconsacc = protacc
                leftcons = None
                leftconsacc = None
                startpos = '1'
                endpos = consrowend
                break
            elif j == len(contiglist)-1:
                if consrow < fragrow: # if the last conserved gene on the contig is before the fragment: 
                    leftcons = protname
                    leftconsacc = protacc
                    startpos = consrowstart
                    endpos = str(contig_lengths[contig])
                    rightcons = None
                    rightconsacc = None
                    break
            else: 
                nextconsrow = contiglist[j+1][0]
                nextconsrowstart = gff1[nextconsrow][3]
                nextconsrowend = gff1[nextconsrow][4]
                nextprotname = contiglist[j+1][1]
                nextprotacc = contiglist[j+1][2]
                if consrow < fragrow < nextconsrow: 
                    leftcons = protname
                    leftconsacc = protacc
                    rightcons = nextprotname
                    rightconsacc = nextprotacc
                    startpos = consrowstart
                    endpos = nextconsrowend
                    break
        if startpos == None and endpos == None: 
            startpos = '1'
            endpos = endpos = str(contig_lengths[contig])            
        # for every fragment, write the interval coordinates to an output file
        # that can be used to extract the coordinate with esl-sfetch. 
        # the name of the sequence is the fragment name and ITS coordinates (just the fragment). 
        outfile1.write(eslid + '\t' + startpos + '\t' + endpos + '\t' + contig + '\n') 
        print('wrote native')
        # first, for fragments that have two flanking genes (meaning there's a strong way to find them in an outgroup), 
        # identify the corresponding syntenic region in the outgroup.
        # if there's an annotated fragment there with the same name as the present one, call it conserved
        print(leftcons, rightcons)
        if leftcons != None and rightcons != None: 
            focal_flanks = True
            for contig2 in gff2_contig_main_dict: 
                if target_flanks_syntenic or synteny_found: break
                contiglist = [entry for entry in gff2_contig_main_dict[contig2] if entry[1] != fragid] # used to be contiglist = gff2_contig_main_dict[contig2] - this saves against fragments annotated as conserved in other gff by removing them from interfering with finding the real syntenic region
                for j in range(0, len(contiglist)-1): 
                    if target_flanks_syntenic or synteny_found: break
                    consrow = contiglist[j][0]
                    startpos = gff2[consrow][3] # start position is start of leftmost cons 
                    protname = contiglist[j][1]
                    protacc = contiglist[j][2]
                    nextconsrow = contiglist[j+1][0]
                    endpos = gff2[nextconsrow][4] # end position is end of rightmost cons 
                    nextprotname = contiglist[j+1][1]
                    nextprotacc = contiglist[j+1][2]
                    if protname == leftcons and leftcons != rightcons: # if flanking genes are different, ie not in intron
                        if nextprotname == rightcons: 
                            target_flanks_syntenic = True
                            int_rows = gff2[consrow+1:nextconsrow]
                            for int_row in int_rows:
                                if any(s in int_row[8] for s in frag_colors) or any(s in int_row[8] for s in main_colors):
                                    newfragname = extract_name(int_row[8])
                                    newfragacc = extract_acc(int_row[8])
                                    incluster = False
                                    if newfragacc in cluster_dict.get(frag_acc, []):
                                        incluster = True
                                        print(newfragacc, ', in an intermediate row between leftcons/rightcons, in cluster or same - synteny found') 
                                    if newfragname == fragid or incluster == True:
                                        synteny_found = True
                                        break
                            if synteny_found == False:
                                outfile2.write(eslid + '\t' + startpos + '\t' + endpos + '\t' + contig2 + '\n') 
                                print('wrote target')
                                break
                        elif nextprotacc in cluster_dict.get(frag_acc, []) or nextprotname == fragid:
                            print(nextprotacc, ' "cons" flanking gene in cluster with fragment or same as - synteny true')
                            synteny_found = True
                            break
                    elif protname == rightcons and leftcons != rightcons: # if flanking genes are different, ie not in intron
                        if nextprotname == leftcons: 
                            target_flanks_syntenic = True
                            int_rows = gff2[consrow+1:nextconsrow]
                            for int_row in int_rows:
                                if any(s in int_row[8] for s in frag_colors) or any(s in int_row[8] for s in main_colors):
                                    newfragacc = extract_acc(int_row[8])
                                    newfragname = extract_name(int_row[8])
                                    incluster = False
                                    if newfragacc in cluster_dict.get(frag_acc, []):
                                        incluster = True
                                        print(newfragacc, ', in an intermediate row between leftcons/rightcons, in cluster or same - synteny found') 
                                    if newfragname == fragid or incluster == True:
                                        synteny_found = True
                                        break
                            if synteny_found == False:
                                outfile2.write(eslid + '\t' + startpos + '\t' + endpos + '\t' + contig2 + '\n') 
                                print('wrote target')
                                break
                        elif nextprotacc in cluster_dict.get(frag_acc, [])  or nextprotname == fragid:
                            print(nextprotacc, ' "cons" flanking gene in cluster with fragment or same as - synteny true')
                            synteny_found = True
                            break
                    elif protname == rightcons == leftcons: # if intron, take the whole range of the parent gene; allows for gene structure changes/annotation errors
                        currprotname = rightcons
                        currprotacc = rightconsacc
                        added = 1
                        while currprotname == rightcons and j+1+added < len(contiglist):
                            currprotname = contiglist[j+1+added][1]
                            currprotacc = contiglist[j+1+added][2]
                            nextconsrow = contiglist[j+1+added][0]
                            endpos = gff2[nextconsrow][4]
                            added += 1
                            int_rows = gff2[consrow+1:nextconsrow]
                            for int_row in int_rows:
                                if any(s in int_row[8] for s in frag_colors) or any(s in int_row[8] for s in main_colors): # these guard against differential main/frag annotaiton in different gffs
                                    newfragacc = extract_acc(int_row[8])
                                    newfragname = extract_name(int_row[8])
                                    incluster = False
                                    if newfragacc in cluster_dict.get(frag_acc, []):
                                        incluster = True
                                        print(newfragacc, ', intron case, in an intermediate row between leftcons/rightcons, in cluster or same - synteny found') 
                                    if newfragname == fragid or incluster == True:
                                        synteny_found = True
                                        break
                        if currprotname != rightcons and synteny_found == False:
                            outfile2.write(eslid + '\t' + startpos + '\t' + endpos + '\t' + contig2 + '\n')  
                            print('wrote target')
                            target_flanks_syntenic = True
                            break
                        elif currprotacc in cluster_dict.get(frag_acc, []):
                            print(currprotacc, ' intron case,a in cluster with fragment or same as - synteny true')
                            synteny_found = True
                            target_flanks_syntenic = True
                            break
        if (leftcons != None and rightcons == None) or (leftcons != None and rightcons != None and synteny_found == False): # could be because target_flanks_syntenic = False or because target_flanks_syntenic = True and just not in region
            # if either flank is missing or the between-flank search failed, search on either side of each flank
            # start with leftcons 
            leftcons_found = False
            for contig2 in gff2_contig_main_dict: 
                if leftcons_found or synteny_found: break # stop iterating over contigs once you've found what you're looking for
                contiglist = [entry for entry in gff2_contig_main_dict[contig2] if entry[1] != fragid] # used to be contiglist = gff2_contig_main_dict[contig2] - this saves against fragments annotated as conserved in other gff by removing them from interfering with finding the real syntenic region
                for j in range(0, len(contiglist)): # used t be -1 
                    if leftcons_found or synteny_found: break                        
                    consrow = contiglist[j][0]
                    #startpos = gff2[consrow][3] # start position is start of leftmost cons 
                    protname = contiglist[j][1]
                    #endpos = gff2[nextconsrow][4] # end position is end of rightmost cons 
                    if protname == leftcons: 
                        # right bound
                        rightidx = 0
                        newgeneright = False
                        while newgeneright == False and j+1+rightidx < len(contiglist):
                            nextprotname = contiglist[j+1+rightidx][1]
                            if nextprotname != protname: 
                                newgeneright = True
                            else: 
                                rightidx += 1
                        if newgeneright == False: 
                            nextconsrow = consrow + 1
                            while nextconsrow < len(gff2) and gff2[nextconsrow][0] == contig2:
                                nextconsrow += 1 
                            endpos = gff2[nextconsrow - 1][4] # nextconsrow, endpos is now the first row index NOT on contig2 (or len(gff2))
                        else:
                            nextconsrow = contiglist[j+1+rightidx][0]
                            endpos = gff2[nextconsrow][4] # end position is end of rightmost cons 
                        # left bound
                        if j == 0:
                            startpos = '1'
                            prevconsrow = consrow - 1
                        else:
                            leftidx = 0
                            newgeneleft = False
                            while newgeneleft == False and j-1-leftidx >= 0:
                                lastprotname = contiglist[j-1-leftidx][1]
                                if lastprotname != protname: 
                                    newgeneleft = True
                                else: 
                                    leftidx += 1
                            if newgeneleft == False:
                                prevconsrow = consrow - 1
                                while prevconsrow >= 0 and gff2[prevconsrow][0] == contig2:
                                    prevconsrow -= 1 # prevconsrow is now the last row index NOT on contig2 (or -1)
                                startpos = gff2[prevconsrow + 1][3]
                            else:
                                prevconsrow = contiglist[j-1-leftidx][0]
                                startpos = gff2[prevconsrow][3] # start position is start of leftmost cons 
                        int_rows = gff2[prevconsrow+1:nextconsrow] # take conserved genes on either side of single flank and use them as boundaries
                        for int_row in int_rows:
                            if any(s in int_row[8] for s in frag_colors) or any(s in int_row[8] for s in main_colors):
                                newfragacc = extract_acc(int_row[8])
                                newfragname = extract_name(int_row[8])
                                incluster = False
                                if newfragacc in cluster_dict.get(frag_acc, []):
                                    incluster = True
                                if newfragname == fragid or incluster == True:
                                    synteny_found = True
                                    break
                        if synteny_found == False and (focal_flanks==False or target_flanks_syntenic == False): # if only one flank OR synteny broken in outgroup
                            outfile2.write(eslid + '-L' + '\t' + startpos + '\t' + endpos + '\t' + contig2 + '\n') 
                        leftcons_found = True
                        break
        if (leftcons == None and rightcons != None) or (leftcons != None and rightcons != None and synteny_found == False):
            rightcons_found = False
            for contig2 in gff2_contig_main_dict: 
                if rightcons_found or synteny_found: break
                contiglist = [entry for entry in gff2_contig_main_dict[contig2] if entry[1] != fragid] # used to be contiglist = gff2_contig_main_dict[contig2] - this saves against fragments annotated as conserved in other gff by removing them from interfering with finding the real syntenic region
                for j in range(0, len(contiglist)): # used to be -1
                    if rightcons_found or synteny_found: break                        
                    consrow = contiglist[j][0]
                    #startpos = gff2[consrow][3] # start position is start of leftmost cons 
                    protname = contiglist[j][1]
                    #endpos = gff2[nextconsrow][4] # end position is end of rightmost cons 
                    if protname == rightcons: 
                        # right bound
                        rightidx = 0
                        newgeneright = False
                        while newgeneright == False and j+1+rightidx < len(contiglist):
                            nextprotname = contiglist[j+1+rightidx][1]
                            if nextprotname != protname: 
                                newgeneright = True
                            else: 
                                rightidx += 1
                        if newgeneright == False: 
                            nextconsrow = consrow + 1
                            while nextconsrow < len(gff2) and gff2[nextconsrow][0] == contig2:
                                nextconsrow += 1 
                            endpos = gff2[nextconsrow - 1][4] # nextconsrow, endpos is now the first row index NOT on contig2 (or len(gff2))
                        else:
                            nextconsrow = contiglist[j+1+rightidx][0]
                            endpos = gff2[nextconsrow][4] # end position is end of rightmost cons 
                        # left bound
                        if j == 0:
                            startpos = '1'
                            prevconsrow = consrow - 1
                        else:
                            leftidx = 0
                            newgeneleft = False
                            while newgeneleft == False and j-1-leftidx >= 0:
                                lastprotname = contiglist[j-1-leftidx][1]
                                if lastprotname != protname: 
                                    newgeneleft = True
                                else: 
                                    leftidx += 1
                            if newgeneleft == False:
                                prevconsrow = consrow - 1
                                while prevconsrow >= 0 and gff2[prevconsrow][0] == contig2:
                                    prevconsrow -= 1 # prevconsrow is now the last row index NOT on contig2 (or -1)
                                startpos = gff2[prevconsrow + 1][3]
                            else:
                                prevconsrow = contiglist[j-1-leftidx][0]
                                startpos = gff2[prevconsrow][3] # start position is start of leftmost cons 
                        int_rows = gff2[prevconsrow+1:nextconsrow] # take conserved genes on either side of single flank and use them as boundaries
                        for int_row in int_rows:
                            if any(s in int_row[8] for s in frag_colors) or any(s in int_row[8] for s in main_colors):
                                newfragname = extract_name(int_row[8])
                                newfragacc = extract_acc(int_row[8])
                                incluster = False
                                if newfragacc in cluster_dict.get(frag_acc, []):
                                    incluster = True
                                if newfragname == fragid or incluster == True:
                                    synteny_found = True
                                    break
                        if synteny_found == False and (focal_flanks==False or target_flanks_syntenic == False):
                            outfile2.write(eslid + '-R' + '\t' + startpos + '\t' + endpos + '\t' + contig2 + '\n') 
                        rightcons_found = True
                        break
        if focal_flanks == True and target_flanks_syntenic == False: # keep track of fragments with solid synteny in focal that is lost in target
            outfile3.write(eslid + '\n')
        outfile4.write(eslid + '\t' + str(synteny_found) + '\n')

        
outfile1.close()
outfile2.close()
outfile3.close()  
outfile4.close()
            
    
