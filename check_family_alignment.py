#!/usr/bin/env python3
"""
Print a pairwise %identity matrix for one gene family's alignment.
Usage: python3 check_family_alignment.py <Core|Shell> <family_id>

What to look for:
  - Smoothly varying identity across all pairs (e.g. all in the 90-99%
    range) -> looks like real population-level sequence diversity.
  - A sharp two-block pattern (some pairs ~99% identical, others only
    ~70-80%, with a clear split between two subsets of sequences)
    -> classic signature of two paralogous genes accidentally merged
    into one "family" by clustering, not real within-family diversity.
"""

import sys
import os

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
GAP_CHARS = {"-", "."}


def parse_fasta_aln(path):
    seqs = {}
    name = None
    cur = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(cur)
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
        if name:
            seqs[name] = "".join(cur)
    return seqs


def pct_identity(a, b):
    compared = 0
    same = 0
    for x, y in zip(a, b):
        if x in GAP_CHARS or y in GAP_CHARS:
            continue
        compared += 1
        if x.upper() == y.upper():
            same += 1
    return (100 * same / compared) if compared > 0 else float("nan")


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 check_family_alignment.py <Core|Shell> <family_id>")
        sys.exit(1)

    category, family_id = sys.argv[1], sys.argv[2]
    if category not in FOLDERS:
        print(f"category must be 'Core' or 'Shell', got '{category}'")
        sys.exit(1)

    path = os.path.join(FOLDERS[category], family_id + ".aln")
    if not os.path.exists(path):
        print(f"File not found: {path}")
        sys.exit(1)

    seqs = parse_fasta_aln(path)
    names = list(seqs.keys())
    aln_len = len(next(iter(seqs.values())))

    print(f"Family: {family_id} | {len(names)} sequences | alignment length {aln_len}\n")
    print("Pairwise %identity matrix:\n")

    header = "          " + "".join(f"{n[:9]:>10}" for n in names)
    print(header)
    for n1 in names:
        row = f"{n1[:9]:>10}"
        for n2 in names:
            pid = pct_identity(seqs[n1], seqs[n2])
            row += f"{pid:10.1f}"
        print(row)


if __name__ == "__main__":
    main()
