# load a) blast file, b) main hit file, and c) bacterial hit file
# append a column to a) and c), blast and bacterial dfs, to specify source, 'native' or 'bacterial'
# concatenate a) and c) and then sort by contig and then position
# for each hit in concatenated df, check if it overlaps (or is identical to) any of the main hits
# if it is not, check whether there's a hit in the existing hit vector that overlaps
# if there is overlap, determine which to keep, based only on e-value

###### VERSION INFO
## This version doesn't have a separate criterion for bacterial sequences; just uses which e-value is smaller.
## updated 7/23/25 to fix bug in logic to prevent overlapping extra hits 
## updated 7/30/25 to fix bug allowing multiple overlapping extra hits with len < overlap_len 

import sys
import numpy as np 
# set parameters
overlap_limit = 50 # maximum allowed overlap between hits

# load files and do checks 
if len(sys.argv) != 5:
    sys.stderr.write(f"Usage: {sys.argv[0]} <full_blast> <bacterial_blast> <TE blast> <main_hits>\n")
    sys.exit(1) 

full_blast_path = sys.argv[1]
bacterial_blast_path = sys.argv[2]
TE_blast_path = sys.argv[3]
main_hits_path = sys.argv[4]

full_fly_blast = np.genfromtxt(full_blast_path, delimiter="\t", dtype=str)
full_bacterial_blast = np.genfromtxt(bacterial_blast_path, delimiter="\t", dtype=str)
full_TE_blast = np.genfromtxt(TE_blast_path, delimiter="\t", dtype=str)
main_hits = np.genfromtxt(main_hits_path, delimiter="\t", dtype=str)

if full_fly_blast.shape[1] != full_bacterial_blast.shape[1] or full_fly_blast.shape[1] != full_TE_blast.shape[1] or full_fly_blast.shape[1] != main_hits.shape[1]:
    sys.stderr.write("Error: Full fly blast, bacterial blast, TE blast, and main hits files have different number of columns.\n")
    sys.exit(1)

# define output files
out_file = open(full_blast_path + "_extra_hits", 'w')

# specify columns - must be the same for both blast files
COLUMNS = [
    "prot_id", "prot_len", "prot_start", "prot_end",
    "contig_id", "contig_len",
    "contig_start", "contig_end",
    "evalue", "pident", "length", "frame"
    ]
prot_id_col = COLUMNS.index("prot_id")
contig_id_col = COLUMNS.index("contig_id")
contig_start_col = COLUMNS.index("contig_start")
contig_end_col = COLUMNS.index("contig_end")
evalue_col = COLUMNS.index("evalue")

# concatenate full fly and bacterial blast files
concatenated_full_blast = np.concatenate((full_fly_blast, full_bacterial_blast, full_TE_blast), axis=0)
# Sort concatenated full blast
# first by contig id and then by start position (lexsort sorts by the last column first)
sorting_indices = np.lexsort((
    concatenated_full_blast[:, contig_start_col].astype(int),  # secondary sort
    concatenated_full_blast[:, contig_id_col]                  # primary sort
))
concatenated_full_blast = concatenated_full_blast[sorting_indices].tolist() ###!

# sort main hits by contig id and then by start position, same as aboves
sorting_indices_main_hits = np.lexsort((
    main_hits[:, contig_start_col].astype(int),  # secondary sort
    main_hits[:, contig_id_col]                  # primary sort
))
main_hits = main_hits[sorting_indices_main_hits].tolist() ###!

# create dictionary for main hits for quick lookup
main_hits_dict = {}
for hit in main_hits:
    contig_id = hit[contig_id_col]
    if contig_id not in main_hits_dict:
        main_hits_dict[contig_id] = []
    main_hits_dict[contig_id].append(hit)

accepted_new_hits_dict = {}
for idx, row in enumerate(concatenated_full_blast):
    if idx % 10 == 0:
        print(f"Processing row index: {idx} out of {len(concatenated_full_blast)}")
    contig_id = row[contig_id_col]
    contig_start = min(int(row[contig_start_col]), int(row[contig_end_col]))
    contig_end = max(int(row[contig_start_col]), int(row[contig_end_col]))
    hit_len = contig_end - contig_start
    main_overlaps = False
    # Check if this hit overlaps with any main hits
    if contig_id in main_hits_dict:
        for main_hit in main_hits_dict[contig_id]:
            main_start = min(int(main_hit[contig_start_col]), int(main_hit[contig_end_col]))
            main_end = max(int(main_hit[contig_start_col]), int(main_hit[contig_end_col]))
            # Check for overlap
            overlap_start = max(contig_start, main_start)
            overlap_end = min(contig_end, main_end)
            overlap_len = overlap_end - overlap_start
            if overlap_len >= overlap_limit or overlap_len >= hit_len:
                main_overlaps = True
                break
            else: 
                continue
    # if no overlaps with main hits, check for overlaps with accepted new hits
    if main_overlaps == False:
        to_remove = []
        new_accepted_overlaps = False
        if contig_id in accepted_new_hits_dict: 
            for accepted_hit in accepted_new_hits_dict[contig_id][-30:]:
                acc_start = min(int(accepted_hit[contig_start_col]), int(accepted_hit[contig_end_col]))
                acc_end = max(int(accepted_hit[contig_start_col]), int(accepted_hit[contig_end_col]))
                # Check for overlap with accepted hits and assess which hit to keep
                # if no overlap, add to accepted new hits
                overlap_start = max(contig_start, acc_start)
                overlap_end = min(contig_end, acc_end)
                overlap_len = overlap_end - overlap_start
                if overlap_len < overlap_limit and overlap_len < hit_len:
                    continue
                else:
                    if float(row[evalue_col]) < float(accepted_hit[evalue_col]):
                        #accepted_new_hits_dict[contig_id].remove(accepted_hit)
                        #accepted_new_hits_dict[contig_id].append(row)
                        #appended = True
                        to_remove.append(accepted_hit)
                    else: 
                        new_accepted_overlaps = True
                        break
            if new_accepted_overlaps == False:
                # If no overlaps with accepted hits, add to accepted new hits
                for entry in to_remove:
                    accepted_new_hits_dict[contig_id].remove(entry)
                accepted_new_hits_dict[contig_id].append(row)
        else: 
            # If no accepted hits for this contig, create a new entry
            accepted_new_hits_dict[contig_id] = [row]
            
for contig_id, hits in accepted_new_hits_dict.items():
    for hit in hits:
        out_file.write("\t".join(hit) + "\n")

out_file.close()
