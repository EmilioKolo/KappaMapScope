#!/usr/bin/env python3

import argparse
from mappability_core import compute_mappability


def main():
    parser = argparse.ArgumentParser(
        description="Compute mappability for a single k."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-k", "--kmer", type=int, required=True)
    parser.add_argument("-x", "--index", required=True)
    parser.add_argument("--output-prefix", default="mappability")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--cleanup", action="store_true", 
                        help="Remove intermediate SAM files after processing")
    parser.add_argument("--create-bam", action="store_true", 
                        help="Create sorted BAM file from SAM output")

    args = parser.parse_args()

    print("Computing mappability...")
    m = compute_mappability(
        args.input,
        args.kmer,
        args.index,
        args.output_prefix,
        threads=args.threads,
        make_bam=args.create_bam,
        cleanup=args.cleanup
    )

    out_file = f"{args.output_prefix}.k{args.kmer}.txt"

    with open(out_file, "w") as out:
        out.write(f"k\tmappability\n")
        out.write(f"{args.kmer}\t{m:.6f}\n")

    print(f"Mappability: {m:.6f}")
    print(f"Output: {out_file}")


if __name__ == "__main__":
    main()