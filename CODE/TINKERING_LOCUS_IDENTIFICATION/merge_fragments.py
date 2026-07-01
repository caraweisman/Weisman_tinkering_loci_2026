import numpy as np
import sys
import re 
import math
import warnings
warnings.filterwarnings('ignore')

# 6-27-26 version

dist_thresh = 10000
prot_dist_thresh = 200

extra_hits = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', comments=None)

clusters = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t', invalid_raise=False, filling_values='')
contig_col = 0
contig_start_col = 3
contig_end_col = 4
strand_col = 6
desc_col = 8

cluster_dict = {}
with open(sys.argv[2]) as f:
    for line in f:
        fields = line.strip().split('\t')
        focal = fields[0]
        hits = fields[1:]
        cluster_dict[focal] = hits


def parse_description(description):

    # (130-184/292)
    m = re.search(r"ID=\((\d+)-(\d+)/(\d+)\)-", description)
    protstart = int(m.group(1))
    protend   = int(m.group(2))
    protlen   = int(m.group(3))

    # e-value
    m = re.search(r"description=([^;]+)", description)
    e_str = re.sub(r"^evalue_", "", m.group(1))
    evalue = float(e_str)

    # capture protein ID through species tag
    m = re.search(r"\b([NXY]P_\d+(?:\.\d+)?[^;]*?\[Drosophila[- ]melanogaster\])", description)
    protID = m.group(1)

    # remove isoform chunk only
    protID = re.sub(r",?-isoform\b[^;\[]*", "", protID)

    # ensure exactly ONE dash before species bracket
    protID = re.sub(
        r"-*\[(Drosophila[- ]melanogaster)\]",
        r"-[\1]",
        protID
    )

    return protID, float(evalue), int(protlen), int(protstart), int(protend)

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

def extract_accession(desc):
    accession = re.split('-', desc)[2]
    return accession 

output_rows = []

for row in extra_hits: 
    contig = row[contig_col]
    contig_start = min(int(row[contig_start_col]), int(row[contig_end_col]))
    contig_end = max(int(row[contig_start_col]), int(row[contig_end_col]))
    strand = row[strand_col]
    desc = row[desc_col]
    isoformfree_protid, evalue, prot_len, prot_start, prot_end = parse_description(desc)
    name = extract_name(desc)
    #print(name)
    acc = extract_accession(desc)
    if len(output_rows) == 0: 
        output_rows.append(row)
    else: 
        old_contig = output_rows[-1][contig_col]
        old_contig_start = min(int(output_rows[-1][contig_start_col]), int(output_rows[-1][contig_end_col]))
        old_contig_end = max(int(output_rows[-1][contig_start_col]), int(output_rows[-1][contig_end_col]))
        old_strand = output_rows[-1][strand_col]
        old_desc = output_rows[-1][desc_col]
        old_name = extract_name(old_desc)
        old_acc = extract_accession(old_desc)
        if acc in cluster_dict[old_acc] or old_acc in cluster_dict[acc]: 
            in_cluster = True
        else: 
            in_cluster = False
        isoformfree_oldid, old_evalue, old_prot_len, old_prot_start, old_prot_end = parse_description(old_desc)
        correct_order = False
        if old_strand == strand == '+': 
            if prot_start > old_prot_end or prot_start > old_prot_end-prot_dist_thresh:
                correct_order = True
        if old_strand == strand == '-':
            if prot_start  < old_prot_end or prot_start - prot_dist_thresh < old_prot_end:
                correct_order = True
		#prot_gap = min(abs(prot_start - old_prot_start), abs(prot_start - old_prot_end), abs(prot_end - old_prot_start), abs(prot_end - old_prot_end))
        if (name == old_name or in_cluster == True) and contig_start - old_contig_end < dist_thresh and contig == old_contig and strand == old_strand and correct_order == True:
            print(row)
            print(output_rows[-1])
            combined_row = []
            combined_row.append(contig)
            combined_row.append('tblastn') 
            combined_row.append('match') 
            newstart = min(old_contig_start, contig_start)
            combined_row.append(newstart)
            newend = max(old_contig_end, contig_end)
            combined_row.append(newend)
            combined_row.append('.') 
            combined_row.append(strand)
            combined_row.append('.')
            protstart = min(prot_start, old_prot_start)
            protend = max(prot_end, old_prot_end) 
            color_match = re.search(r'(color=[^;]+)', old_desc)
            namestr = '(' + str(protstart) + '-' + str(protend) + '/' + str(prot_len) + ')-' + isoformfree_protid + '_' + str(newstart) + '_' + str(newend)
            desc = 'ID=' + namestr + ';' + 'Name=' + namestr + ';description=' + str(old_evalue) + ';' + color_match.group() if color_match else ''
            combined_row.append(desc)
            output_rows.pop()
            output_rows.append(combined_row)
        else: 
            output_rows.append(row)

outfile = open(sys.argv[1] + '_MERGEDCLUST_10KB', 'w')
for row in output_rows: 
    outfile.write('\t'.join(map(str, row)))
    outfile.write('\n')
outfile.close()


