#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from mappability_core import compute_mappability
from multiprocessing import Pool


def compute_wrapper(args_tuple):
    """
    Wrapper to unpack arguments for multiprocessing.
    """
    fasta, k, index, prefix, threads = args_tuple
    m = compute_mappability(
        fasta, k, index, prefix,
        threads=threads, cleanup=True
    )
    return (k, m)


def parse_k_values(k_list, k_range):
    if k_list:
        return [int(k) for k in k_list.split(",")]
    elif k_range:
        start, end, step = map(int, k_range.split(","))
        return list(range(start, end + 1, step))
    else:
        raise ValueError("Provide --k-list or --k-range")


def split_fasta(input_fasta, out_prefix):
    """
    Split multi-FASTA into one FASTA per chromosome.
    Returns dict: {chrom_name: fasta_path}
    """
    chrom_files = {}

    current_name = None
    current_lines = []

    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_name:
                    path = f"{out_prefix}.{current_name}.fa"
                    with open(path, "w") as out:
                        out.write(f">{current_name}\n")
                        out.writelines(current_lines)
                    chrom_files[current_name] = path

                current_name = line.strip()[1:].split()[0]
                current_lines = []
            else:
                current_lines.append(line)

        # last chromosome
        if current_name:
            path = f"{out_prefix}.{current_name}.fa"
            with open(path, "w") as out:
                out.write(f">{current_name}\n")
                out.writelines(current_lines)
            chrom_files[current_name] = path

    return chrom_files


def main():
    parser = argparse.ArgumentParser(
        description="Compute mappability per chromosome over multiple k values."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-x", "--index", required=True)
    parser.add_argument("--k-list",
                        help="Comma-separated k values (e.g. 21,31,41)")
    parser.add_argument("--k-range",
                        help="start,end,step (e.g. 20,100,10)")
    parser.add_argument("--output-prefix", default="mappability_sweep")
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="Number of parallel processes")
    parser.add_argument("-t", "--threads", type=int, default=1)

    args = parser.parse_args()

    ks = parse_k_values(args.k_list, args.k_range)

    print("Splitting FASTA by chromosome...")
    chrom_fastas = split_fasta(args.input, args.output_prefix)

    print(f"Found {len(chrom_fastas)} chromosomes.")

    for chrom, fasta_path in chrom_fastas.items():
        print(f"\n[Chromosome: {chrom}] Running k sweep...")

        prefix = f"{args.output_prefix}.{chrom}"

        tasks = [
            (fasta_path, k, args.index, prefix, args.threads)
            for k in ks
        ]

        if args.processes > 1:
            with Pool(args.processes) as pool:
                results = pool.map(compute_wrapper, tasks)
        else:
            results = [compute_wrapper(t) for t in tasks]

        results.sort(key=lambda x: x[0])

        # Write table
        table_file = f"{prefix}.table.tsv"
        with open(table_file, "w") as out:
            out.write("k\tmappability\n")
            for k, m in results:
                out.write(f"{k}\t{m:.6f}\n")

        # Plot
        ks_vals = [k for k, _ in results]
        m_vals = [m for _, m in results]

        plt.figure()
        plt.plot(ks_vals, m_vals, marker='o')
        plt.xlabel("k-mer size")
        plt.ylabel("Mappability")
        plt.title(f"Mappability vs k ({chrom})")
        plt.grid(True)

        fig_file = f"{prefix}.png"
        plt.savefig(fig_file)

        print(f"[{chrom}] Done.")
        print(f"Table: {table_file}")
        print(f"Figure: {fig_file}")

    print("\nAll chromosomes processed.")


if __name__ == "__main__":
    main()