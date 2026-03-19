#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt


def read_table(path):
    """
    Read k vs mappability table.
    """
    ks = []
    ms = []

    with open(path) as f:
        next(f)  # skip header
        for line in f:
            k, m = line.strip().split("\t")
            ks.append(int(k))
            ms.append(float(m))

    return ks, ms


def main():
    parser = argparse.ArgumentParser(description="Plot mappability vs k from table.")
    parser.add_argument("-i", "--input", required=True, help="Input table (TSV)")
    parser.add_argument("-o", "--output", default="mappability_plot.png")
    parser.add_argument("--xlim", help="xmin,xmax")
    parser.add_argument("--ylim", help="ymin,ymax")
    parser.add_argument("--logx", action="store_true", help="Use log scale on x axis")

    args = parser.parse_args()

    ks, ms = read_table(args.input)

    plt.figure()
    plt.plot(ks, ms, marker='o')
    plt.xlabel("k-mer size")
    plt.ylabel("Mappability")
    plt.title("Mappability vs k")
    plt.grid(True)

    # ---- axis controls ----
    if args.xlim:
        xmin, xmax = map(float, args.xlim.split(","))
        plt.xlim(xmin, xmax)

    if args.ylim:
        ymin, ymax = map(float, args.ylim.split(","))
        plt.ylim(ymin, ymax)

    if args.logx:
        plt.xscale("log")

    plt.savefig(args.output)

    print(f"Plot saved to: {args.output}")


if __name__ == "__main__":
    main()