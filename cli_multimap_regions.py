#!/usr/bin/env python3

import argparse
import pysam


def merge_intervals(intervals):
    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or merged[-1][1] < start:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Extract multi-mapping regions from BAM"
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM")
    parser.add_argument("-k", "--kmer", type=int, required=True)
    parser.add_argument("-o", "--output", default="multimapping_regions.bed")

    args = parser.parse_args()

    print("[STEP 1] Counting alignments per k-mer...")
    counts = {}

    with pysam.AlignmentFile(args.input, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            q:str = str(read.query_name)
            counts[q] = counts.get(q, 0) + 1

    print("[STEP 2] Collecting multi-mapping intervals...")
    intervals = []

    with pysam.AlignmentFile(args.input, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            q:str = str(read.query_name)
            if counts.get(q, 0) <= 1:
                continue  # skip unique

            # extract k-mer index
            try:
                idx = int(q.split("_")[1])
            except:
                continue

            start = idx
            end = idx + args.kmer

            intervals.append((start, end))

    print("[STEP 3] Merging intervals...")
    merged = merge_intervals(intervals)

    print("[STEP 4] Writing BED...")
    with open(args.output, "w") as out:
        for start, end in merged:
            out.write(f"contig\t{start}\t{end}\n")

    print(f"Done. Regions written to: {args.output}")
    print(f"Total regions: {len(merged)}")

if __name__ == "__main__":
    main()