#!/usr/bin/env python3
"""
Count specific feature types in the GFF3 files for LyasA and LyasB.
Prints every unique feature type actually present (diagnostic, since
exact Bakta naming/capitalization for some less-common types isn't
being assumed), then counts the specifically requested set.
"""

import os
from collections import Counter

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
GENOMES = ["LyasA", "LyasB"]

REQUESTED_TYPES = [
    "CDS", "tRNA", "tmRNA", "rRNA", "ncRNA", "ncRNA_region",
    "CRISPR", "sORF", "gap", "oriC", "oriV", "oriT",
]


def count_features(path):
    counter = Counter()
    with open(path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            counter[cols[2]] += 1
    return counter


def main():
    all_counts = {}
    for genome in GENOMES:
        path = os.path.join(GFF3_DIR, f"{genome}.gff3")
        counts = count_features(path)
        all_counts[genome] = counts

        print(f"=== {genome}: all unique feature types found ===")
        for ftype, n in sorted(counts.items(), key=lambda x: -x[1]):
            print(f"  {ftype}: {n}")
        print()

    print("=== Requested feature type counts ===")
    header = f"{'Feature':<15}" + "".join(f"{g:>10}" for g in GENOMES)
    print(header)
    print("-" * len(header))
    for ftype in REQUESTED_TYPES:
        row = f"{ftype:<15}"
        for genome in GENOMES:
            # case-insensitive match against actual types found, in case
            # exact casing differs from the requested spelling
            match_count = 0
            for actual_type, n in all_counts[genome].items():
                if actual_type.lower().replace("_", " ").replace("-", " ") == \
                   ftype.lower().replace("_", " ").replace("-", " "):
                    match_count += n
            row += f"{match_count:>10}"
        print(row)


if __name__ == "__main__":
    main()
