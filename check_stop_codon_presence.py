#!/usr/bin/env python3
"""
Check whether Core family alignment sequences end with a stop codon
(TAA/TAG/TGA) or not -- determines whether stop-gain/stop-loss
classification is even possible with this data source.
"""

import glob
import os
import random

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDER = os.path.join(BASE_DIR, "msa_core", "msa_core_dna")
GAP_CHARS = {"-", "."}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def parse_fasta_aln(path):
    seqs = []
    cur = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append("".join(cur))
    return seqs


def main():
    random.seed(42)
    files = sorted(glob.glob(os.path.join(FOLDER, "*.aln")))
    sample = random.sample(files, min(100, len(files)))

    n_checked = 0
    n_ends_in_stop = 0

    for path in sample:
        seqs = parse_fasta_aln(path)
        for s in seqs:
            ungapped = s.upper().replace("-", "").replace(".", "")
            if len(ungapped) < 3:
                continue
            last_codon = ungapped[-3:]
            n_checked += 1
            if last_codon in STOP_CODONS:
                n_ends_in_stop += 1

    print(f"Checked {n_checked} sequences across {len(sample)} sampled Core families")
    print(f"Sequences ending in a stop codon (TAA/TAG/TGA): {n_ends_in_stop} "
          f"({100*n_ends_in_stop/n_checked:.1f}%)")

    # Show a concrete example either way
    example_path = sample[0]
    example_seqs = parse_fasta_aln(example_path)
    print(f"\nExample from {os.path.basename(example_path)}:")
    for s in example_seqs[:3]:
        ungapped = s.upper().replace("-", "").replace(".", "")
        print(f"  length={len(ungapped)}, last 9bp = ...{ungapped[-9:]}")


if __name__ == "__main__":
    main()
