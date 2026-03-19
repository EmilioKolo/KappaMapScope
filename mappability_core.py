#!/usr/bin/env python3

import os
import pysam
import subprocess


def compute_mappability(
    fasta,
    k,
    index,
    prefix,
    threads=1,
    make_bam=False,
    cleanup=False
):
    """
    Full pipeline for one k.
    """
    seq = read_fasta(fasta)

    kmers_fasta = f"{prefix}.k{k}.kmers.fa"
    sam_file = f"{prefix}.k{k}.sam"

    # k-mer caching
    if os.path.exists(kmers_fasta):
        print(f"[CACHE] Reusing existing k-mer file: {kmers_fasta}")
    else:
        print(f"[BUILD] Generating k-mers: {kmers_fasta}")
        generate_kmers(seq, k, kmers_fasta)

    # Alignment caching
    if os.path.exists(sam_file):
        print(f"[CACHE] Reusing existing alignment: {sam_file}")
    else:
        print(f"[ALIGN] Running Bowtie2 for k={k}")
        run_bowtie2(index, kmers_fasta, sam_file, threads)

    # Create BAM file (if requested)
    bam_file = f"{prefix}.k{k}.sorted.bam"
    if make_bam:
        print(f"[BAM] Converting SAM to sorted BAM: {bam_file}")
        pysam.sort("-o", bam_file, sam_file)
        pysam.index(bam_file)

    unique, total = count_unique_mappings(sam_file)

    # Cleanup intermediate files
    if cleanup:
        if os.path.exists(sam_file):
            print(f"[CLEANUP] Removing SAM file: {sam_file}")
            os.remove(sam_file)

    return unique / total if total else 0.0


def count_unique_mappings(sam_file):
    counts = {}

    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            qname = read.query_name
            counts[qname] = counts.get(qname, 0) + 1

    unique = sum(1 for v in counts.values() if v == 1)
    total = len(counts)

    return unique, total


def generate_kmers(seq, k, out_fasta):
    with open(out_fasta, "w") as out:
        for i in range(len(seq) - k + 1):
            out.write(f">kmer_{i}\n{seq[i:i+k]}\n")


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def run_bowtie2(index, kmers_fasta, sam_out, threads=1):
    cmd = [
        "bowtie2",
        "-x", index,
        "-f", kmers_fasta,
        "-S", sam_out,
        "--very-sensitive",
        "-k", "5",          # Number of accepted multi-alignments
        "-p", str(threads)
    ]
    subprocess.run(cmd, check=True)
