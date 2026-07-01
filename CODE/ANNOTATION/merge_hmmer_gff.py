#!/usr/bin/env python3

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
outfile = sys.argv[3]

def parse(line):
    fields = line.strip().split("\t")
    contig = fields[0]
    start  = int(fields[3])
    end    = int(fields[4])
    if start > end:
        start, end = end, start

    attr = fields[8]
    desc = attr.split("description=")[1].split(";")[0]
    if desc.startswith("evalue_"):
        desc = desc[len("evalue_"):]
    evalue = float(desc)

    return [contig, start, end, evalue, line]


# load into memory
f1 = [parse(l) for l in open(file1)]
f2 = [parse(l) for l in open(file2)]

# sort by contig then start
f1.sort(key=lambda x: (x[0], x[1]))
f2.sort(key=lambda x: (x[0], x[1]))

i = 0
j = 0

while i < len(f1) and j < len(f2):
    contig1, s1, e1, ev1, line1 = f1[i]
    contig2, s2, e2, ev2, line2 = f2[j]

    # different contigs → advance the smaller one
    if contig1 < contig2:
        i += 1
        continue
    if contig1 > contig2:
        j += 1
        continue

    # same contig
    overlap = min(e1, e2) - max(s1, s2) + 1

    if overlap <= 10:
        # not a conflict (including no overlap)
        # advance whichever ends first
        if e1 < e2:
            i += 1
        else:
            j += 1
        continue

    # conflict (>10 bp overlap)
    if ev1 <= ev2:
        f2[j][4] = None
        j += 1
    else:
        f1[i][4] = None
        i += 1


with open(outfile, "w") as out:
    for entry in f1 + f2:
        if entry[4] is not None:
            out.write(entry[4])

