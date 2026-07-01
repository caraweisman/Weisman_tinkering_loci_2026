######### This is version 2, different than the files in the V1 subdirectory; rewritten, but takes generally the same approach.


import sys
import pandas as pd

# adjustable parameters
MAX_GAP = 50000
OVERLAP_LIMIT = 40
MINIMUM_COVERAGE = 0.60
PROGRESS_EVERY = 100

COLUMNS = [
    "prot_id", "prot_len", "prot_start", "prot_end",
    "contig_id", "contig_len",
    "contig_start", "contig_end",
    "evalue", "pident", "length", "frame"
    ]

if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <blast.tsv> <protein_ids.txt>\n")
        sys.exit(1)

blast_path = sys.argv[1]
ids_path = sys.argv[2]

out_file = open(blast_path + "_main_hits", 'w')
out_info_file = open(blast_path + "_cov_stats_main_hits", 'w')

#raw_lines = open(blast_path, "r").read().splitlines()
#with open(blast_path, "rb") as f:
#    raw_lines = f.readlines()   # each item en

df = pd.read_csv(blast_path, sep="\t", names=COLUMNS, dtype='str')
#numeric versions for sorting only
df["evalue_num"] = df["evalue"].astype(float)
df["pident_num"] = df["pident"].astype(float)

with open(ids_path) as f:
    protein_ids = f.read().splitlines()

for i, pid in enumerate(protein_ids):
    print(i, 'out of', len(protein_ids))
    query_rows = df[df["prot_id"] == pid].copy()
    
    if len(query_rows) == 0:
        out_info_file.write(f"{pid}\t0\t0\n")
        continue
    prot_len = int(query_rows.iloc[0]["prot_len"])
    query_rows = query_rows.sort_values(["evalue_num", "pident_num"], ascending=[True, False]).reset_index(drop=False)
    final_cluster = None
    final_rows = None
    fallback_cluster = None
    fallback_rows = None

    for anchor_idx in range(len(query_rows)):
        anchor = query_rows.iloc[anchor_idx]
        frame_val = anchor["frame"].split("/")[-1]
        strand = "+" if int(frame_val) > 0 else "-"
        cluster_rows = [anchor_idx]
        cluster_info = [anchor]
        changed = True
        while changed:
            changed = False
            for j in range(len(query_rows)):
                if j in cluster_rows:
                    continue
                cand_row = query_rows.iloc[j]
                # same contig
                if cand_row["contig_id"] != anchor["contig_id"]:
                    continue
                # same strand
                cand_frame_val = cand_row["frame"].split("/")[-1] 
                cand_strand = "+" if int(cand_frame_val) > 0 else "-"
                if cand_strand != strand:
                    continue
                # within MAX_GAP to any cluster hit
                c_start = int(cand_row["contig_start"])
                c_end   = int(cand_row["contig_end"])
                close = False
                for hit in cluster_info:
                    h_start = int(hit["contig_start"])
                    h_end   = int(hit["contig_end"])
                    for p1 in (c_start, c_end):
                        for p2 in (h_start, h_end):
                            if abs(p1 - p2) <= MAX_GAP:
                                    close = True
                if not close:
                        continue
                # proteins overlap by no more than OVERLAP_LIMIT
                c_pstart = int(cand_row["prot_start"])
                c_pend   = int(cand_row["prot_end"])
                overlap_ok = True
                for hit in cluster_info:
                    h_pstart = int(hit["prot_start"])
                    h_pend = int(hit["prot_end"])
                    overlap = max(0, min(c_pend, h_pend) - max(c_pstart, h_pstart) + 1)
                    if overlap > OVERLAP_LIMIT:
                        overlap_ok = False
                if not overlap_ok:
                    continue
                # directional compatibility on contig
                compatible = True
                for hit in cluster_info:
                    h_pstart = int(hit["prot_start"])
                    h_pend   = int(hit["prot_end"])
                    h_gpos = max(int(hit["contig_start"]), int(hit["contig_end"])) # only use ends
                    c_gpos = max(c_start, c_end)
                    if strand == "+":
                        if c_pend > h_pend and c_gpos < h_gpos:
                            compatible = False
                        elif c_pstart < h_pstart and c_gpos > h_gpos:
                            compatible = False
                    else:
                        if c_pend > h_pend and c_gpos > h_gpos:
                            compatible = False
                        elif c_pstart < h_pstart and c_gpos < h_gpos:
                            compatible = False
                if not compatible:
                    continue

                #add to cluster
                cluster_rows.append(j)
                cluster_info.append(cand_row)
                changed = True
        
        # calculate coverage
        intervals = [(int(hit["prot_start"]), int(hit["prot_end"])) for hit in cluster_info]
        intervals.sort()
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1]:
                merged.append((s, e))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        covered = sum(e - s + 1 for s, e in merged)
        coverage = covered / prot_len

        # if the first cluster has more than min coverage, keep it
        # otherwise, look for others that meet this criterion, and bookmark the first one in case none do
        if coverage >= MINIMUM_COVERAGE:
            final_cluster = cluster_info
            final_rows = cluster_rows
            break
        elif anchor_idx == 0:
            fallback_cluster = cluster_info
            fallback_rows = cluster_rows

    # once all anchor points are tried for a gene, pick the best cluster: either the first one if none are >60% id, or the first one if >60% id, or the first with >60% id
    if final_cluster is None:
            final_cluster = fallback_cluster
            final_rows = fallback_rows

    #recompute coverage for the final cluster
    intervals = [(int(hit["prot_start"]), int(hit["prot_end"])) for hit in final_cluster]
    intervals.sort()
    merged = []
    for s, e in intervals:
        if not merged or s > merged[-1][1]:
            merged.append((s, e))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    covered = sum(e - s + 1 for s, e in merged)
    coverage = covered / prot_len
    
    # compute redundant coverage
    intervals = [(int(hit["prot_start"]), int(hit["prot_end"])) for hit in final_cluster]
    redundant_covered = sum(e - s + 1 for s, e in intervals)
    redundant_ratio = redundant_covered / prot_len


    # write final cluster rows to output
    sorted_cluster = sorted(final_cluster, key=lambda x: int(x['prot_start']))

    for hit in sorted_cluster:
        # Convert row to tab-separated string
        row_str = "\t".join(str(hit[col]) for col in COLUMNS)
        out_file.write(f"{row_str}\n")
        #out_file.write(raw_lines[int(hit["index"])] + "\n")
        #out_file.write(raw_lines[int(hit["index"])])

    # write coverage stats to second output
    out_info_file.write(f"{pid}\t{coverage:.2f}\t{redundant_ratio:.2f}\n")

out_file.close()
out_info_file.close()
print("Done.")
