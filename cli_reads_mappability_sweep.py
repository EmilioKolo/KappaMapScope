#!/usr/bin/env python3

import argparse
import gzip
import subprocess
import pysam
import matplotlib.pyplot as plt
from multiprocessing import Pool


def generate_kmers_from_reads(fastq, k, out_fasta):
    """
    Extract k-mers from FASTQ reads.
    """
    # Detect gzip-compressed FASTQ
    if fastq.endswith(".gz"):
        fin_handle = gzip.open(fastq, "rt")
    else:
        fin_handle = open(fastq, "r")

    with fin_handle as fin, open(out_fasta, "w") as out:
        i = 0
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline() # +
            fin.readline() # qual

            for j in range(len(seq) - k + 1):
                out.write(f">read{i}_kmer{j}\n{seq[j:j+k]}\n")
            i += 1


def read_fasta_lengths(fasta):
    """
    Return total genome length from FASTA.
    """
    total = 0
    seq_len = 0

    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                total += seq_len
                seq_len = 0
            else:
                seq_len += len(line.strip())

    total += seq_len
    return total


def run_bowtie2(index, kmers_fasta, sam_out, threads):
    cmd = [
        "bowtie2",
        "-x", index,
        "-f", kmers_fasta,
        "-S", sam_out,
        "--very-sensitive",
        "-k", "2",
        "-p", str(threads)
    ]
    print(f"[ALIGN] Running Bowtie2 on {kmers_fasta}")
    subprocess.run(cmd, check=True)
    print(f"[ALIGN] Finished Bowtie2 → {sam_out}")


def merge_intervals(intervals):
    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or merged[-1][1] < start:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def compute_unique_coverage(sam_file, reference):
    """
    Compute fraction of genome covered by uniquely mapped k-mers.
    """
    unique_reads = {}
    coverage = {}

    print(f"[COVERAGE] Counting alignments in {sam_file}")

    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            qname = read.query_name
            unique_reads[qname] = unique_reads.get(qname, 0) + 1

    # second pass: mark coverage only for unique mappings
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            qname = read.query_name
            if unique_reads.get(qname, 0) != 1:
                continue

            ref = read.reference_name
            start = read.reference_start
            end = read.reference_end

            if ref not in coverage:
                coverage[ref] = []

            coverage[ref].append((start, end))

    covered_bases = 0
    for ref in coverage:
        merged = merge_intervals(coverage[ref])
        for start, end in merged:
            covered_bases += (end - start)

    print(f"[COVERAGE] Covered bases: {covered_bases}")

    # approximate genome length from SAM header
    genome_length = read_fasta_lengths(reference)

    return covered_bases / genome_length if genome_length else 0.0


def compute_wrapper(args_tuple):
    fastq, k, index, prefix, threads, reference = args_tuple

    print(f"[START] k={k}")

    kmers_fasta = f"{prefix}.k{k}.kmers.fa"
    sam_file = f"{prefix}.k{k}.sam"

    print(f"[STREAM] Streaming k-mers directly to Bowtie2 for k={k}")

    bowtie_cmd = [
        "bowtie2",
        "-x", index,
        "-f", "-",                 # read from stdin
        "-S", sam_file,
        "--very-sensitive",
        "-k", "2",
        "-p", str(threads)
    ]

    proc = subprocess.Popen(
        bowtie_cmd,
        stdin=subprocess.PIPE,
        text=True
    )

    seen_kmers = set()

    # stream kmers
    if fastq.endswith(".gz"):
        fin_handle = gzip.open(fastq, "rt")
    else:
        fin_handle = open(fastq, "r")

    with fin_handle as fin:
        i = 0
        for line in fin:
            header = line
            seq = fin.readline().strip()
            fin.readline()
            fin.readline()

            for j in range(len(seq) - k + 1):
                kmer = seq[j:j+k]
                rc = revcomp(kmer)
                canon = min(kmer, rc)

                if canon in seen_kmers:
                    continue

                seen_kmers.add(canon)
                proc.stdin.write(f">kmer_{len(seen_kmers)}\n{canon}\n")

            i += 1

            if i % 100000 == 0:
                print(f"[STREAM] {i} reads processed")

    print(f"[STREAM] Unique k-mers sent: {len(seen_kmers)}")

    proc.stdin.close()
    proc.wait()

    print(f"[STREAM] Finished streaming k={k}")

    m = compute_unique_coverage(sam_file, reference)

    print(f"[DONE] k={k} → mappability={m:.6f}")

    return (k, m)


def parse_k_values(k_list, k_range):
    if k_list:
        return [int(k) for k in k_list.split(",")]
    elif k_range:
        start, end, step = map(int, k_range.split(","))
        return list(range(start, end + 1, step))
    else:
        raise ValueError("Provide --k-list or --k-range")


def main():
    parser = argparse.ArgumentParser(
        description="Compute genome coverage by uniquely mapped read-derived k-mers."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-x", "--index", required=True)
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--k-list")
    parser.add_argument("--k-range")
    parser.add_argument("--output-prefix", default="reads_mappability")
    parser.add_argument("-p", "--processes", type=int, default=1)
    parser.add_argument("-t", "--threads", type=int, default=1)

    args = parser.parse_args()

    ks = parse_k_values(args.k_list, args.k_range)

    tasks = [
        (args.input, k, args.index, args.output_prefix, args.threads, args.reference)
        for k in ks
    ]

    print(f"[INFO] Running ks: {ks}")
    print(f"[INFO] Processes: {args.processes}, Threads per job: {args.threads}")

    if args.processes > 1:
        with Pool(args.processes) as pool:
            results = pool.map(compute_wrapper, tasks)
    else:
        results = [compute_wrapper(t) for t in tasks]

    results.sort(key=lambda x: x[0])

    # write table
    table_file = f"{args.output_prefix}.table.tsv"
    with open(table_file, "w") as out:
        out.write("k\tcoverage_mappability\n")
        for k, m in results:
            out.write(f"{k}\t{m:.6f}\n")

    # plot
    ks_vals = [k for k, _ in results]
    m_vals = [m for _, m in results]

    plt.figure()
    plt.plot(ks_vals, m_vals, marker='o')
    plt.xlabel("k-mer size")
    plt.ylabel("Genome coverage (unique)")
    plt.title("Read-derived mappability vs k")
    plt.grid(True)

    fig_file = f"{args.output_prefix}.png"
    plt.savefig(fig_file)

    print("Done.")
    print(f"Table: {table_file}")
    print(f"Figure: {fig_file}")


if __name__ == "__main__":
    main()