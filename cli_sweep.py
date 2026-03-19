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
    m = compute_mappability(fasta, k, index, prefix, 
                            threads=threads, cleanup=True)
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
        description="Compute mappability over multiple k values."
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

    print("Running k sweep...")

    tasks = [
        (args.input, k, args.index, args.output_prefix, args.threads)
        for k in ks
    ]

    if args.processes > 1:
        with Pool(args.processes) as pool:
            results = pool.map(compute_wrapper, tasks)
    else:
        results = [compute_wrapper(t) for t in tasks]

    results.sort(key=lambda x: x[0])

    # Write table
    table_file = f"{args.output_prefix}.table.tsv"
    with open(table_file, "w") as out:
        out.write("k\tmappability\n")
        for k, m in results:
            out.write(f"{k}\t{m:.6f}\n")

    # Generate plot
    ks_vals = [k for k, _ in results]
    m_vals = [m for _, m in results]

    plt.figure()
    plt.plot(ks_vals, m_vals, marker='o')
    plt.xlabel("k-mer size")
    plt.ylabel("Mappability")
    plt.title("Mappability vs k")
    plt.grid(True)

    fig_file = f"{args.output_prefix}.png"
    plt.savefig(fig_file)

    print("Done.")
    print(f"Table: {table_file}")
    print(f"Figure: {fig_file}")


if __name__ == "__main__":
    main()