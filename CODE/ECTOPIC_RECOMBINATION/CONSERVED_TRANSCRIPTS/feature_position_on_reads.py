#!/usr/bin/env python3
"""
Convert genomic feature coordinates (from a bedtools intersect output of
BAM vs GFF, with -wb -bed -split) into read-relative coordinates.

For each (read, feature) pair in the intersect file, finds where on the
read sequence the feature-aligned portion sits. Introns collapse out
naturally because get_aligned_pairs() only returns aligned bases.

Usage:
    python feature_positions_on_reads.py <intersect.tsv> <input.bam> <output.tsv>
"""

import sys
import pysam
from collections import defaultdict


def main(intersect_path, bam_path, out_path):
    # --- 1. Parse intersect file ---
    # Layout from `bedtools intersect -abam BAM -b GFF -wb -bed -split`:
    #   cols 1-12  : BED12 of the read alignment block
    #   cols 13-21 : GFF feature (seqid, source, type, start, end, score, strand, frame, attributes)
    # 0-indexed: read_name=3, gff_seqid=12, gff_start=15, gff_end=16, gff_attrs=20
    jobs = defaultdict(list)  # contig -> list of (read_name, fstart, fend, feature_attrs)
    n_lines = 0
    with open(intersect_path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 21:
                sys.stderr.write(
                    f"Warning: line has {len(fields)} fields, expected >=21. "
                    f"Check column layout. Skipping.\n"
                )
                continue
            read_name = fields[3]
            contig    = fields[12]
            fstart    = int(fields[15]) - 1   # GFF 1-based -> 0-based half-open
            fend      = int(fields[16])
            attrs     = fields[20]
            jobs[contig].append((read_name, fstart, fend, attrs))
            n_lines += 1

    sys.stderr.write(f"Parsed {n_lines} intersect lines across {len(jobs)} contigs.\n")

    # --- 2. Process BAM contig by contig ---
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_out = 0
    n_missing_read = 0
    n_no_overlap = 0

    with open(out_path, "w") as out:
        out.write("read\tcontig\tfeat_gstart\tfeat_gend\t"
                  "read_start\tread_end\tspan_genomic\tspan_read\tfeature\n")

        for contig, items in jobs.items():
            # Index primary alignments on this contig by read name
            reads_here = {}
            for r in bam.fetch(contig):
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                reads_here[r.query_name] = r

            # Cache aligned_pairs per read (features often share reads)
            pairs_cache = {}

            for read_name, fstart, fend, attrs in items:
                r = reads_here.get(read_name)
                if r is None:
                    n_missing_read += 1
                    continue

                if read_name not in pairs_cache:
                    pairs_cache[read_name] = {
                        ref: q for q, ref in r.get_aligned_pairs()
                        if ref is not None and q is not None
                    }
                pairs = pairs_cache[read_name]

                hits = [pairs[p] for p in range(fstart, fend) if p in pairs]
                if not hits:
                    n_no_overlap += 1
                    continue

                qstart, qend = min(hits), max(hits) + 1
                if r.is_reverse:
                    L = r.query_length
                    qstart, qend = L - qend, L - qstart

                span_g = fend - fstart
                span_r = qend - qstart

                out.write(
                    f"{read_name}\t{contig}\t{fstart}\t{fend}\t"
                    f"{qstart}\t{qend}\t{span_g}\t{span_r}\t{attrs}\n"
                )
                n_out += 1

    sys.stderr.write(
        f"Wrote {n_out} rows.\n"
        f"Skipped {n_missing_read} (read not found in BAM on that contig).\n"
        f"Skipped {n_no_overlap} (feature span entirely intronic on this read).\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
