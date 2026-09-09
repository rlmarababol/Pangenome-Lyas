#!/usr/bin/env python3
"""
Screen every gene family alignment (Core and Shell) for the paralogy /
misclustering signature: some pairs of sequences near-identical, others
well below intraspecific divergence -- i.e. more than one distinct
sequence type merged into a single "family" by clustering.

Flag rule: a family is flagged if its minimum pairwise %identity is
below MIN_ID_FLOOR while its maximum pairwise %identity is at or above
MAX_ID_CEILING. That combination -- some pairs essentially identical,
others clearly not -- is the signature of merged paralogs, not smooth
population-level diversity (which would show a narrower, more
continuous spread of pairwise identities).
"""

import csv
import glob
import os

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
OUT_CSV = os.path.join(BASE_DIR, "family_identity_screen.csv")

GAP_CHARS = {"-", "."}
MIN_ID_FLOOR = 90.0      # some pair below this...
MAX_ID_CEILING = 97.0    # ...while some other pair is at/above this -> flagged


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
    return (100 * same / compared) if compared > 0 else None


def screen_family(path):
    seqs = parse_fasta_aln(path)
    names = list(seqs.keys())
    n = len(names)
    if n < 2:
        return None

    identities = []
    for i in range(n):
        for j in range(i + 1, n):
            pid = pct_identity(seqs[names[i]], seqs[names[j]])
            if pid is not None:
                identities.append(pid)

    if not identities:
        return None

    min_id = min(identities)
    max_id = max(identities)
    mean_id = sum(identities) / len(identities)
    flagged = (min_id < MIN_ID_FLOOR) and (max_id >= MAX_ID_CEILING)

    return {
        "n_seqs": n,
        "min_identity": round(min_id, 2),
        "max_identity": round(max_id, 2),
        "mean_identity": round(mean_id, 2),
        "flagged_misclustered": flagged,
    }


def main():
    rows = []
    for category, folder in FOLDERS.items():
        files = sorted(glob.glob(os.path.join(folder, "*.aln")))
        print(f"{category}: screening {len(files)} families...")
        for idx, path in enumerate(files, 1):
            family_id = os.path.splitext(os.path.basename(path))[0]
            result = screen_family(path)
            if result is None:
                continue
            rows.append({"category": category, "family_id": family_id, **result})
            if idx % 500 == 0:
                print(f"  ...{category}: {idx}/{len(files)}")

    with open(OUT_CSV, "w", newline="") as f:
        fieldnames = ["category", "family_id", "n_seqs", "min_identity",
                      "max_identity", "mean_identity", "flagged_misclustered"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("\n=== Summary: flagged (likely misclustered) families ===")
    for category in FOLDERS:
        cat_rows = [r for r in rows if r["category"] == category]
        n_flagged = sum(1 for r in cat_rows if r["flagged_misclustered"])
        print(f"{category}: {n_flagged} / {len(cat_rows)} families flagged "
              f"({100*n_flagged/len(cat_rows):.1f}%)")

    print(f"\nSaved: {OUT_CSV}")


if __name__ == "__main__":
    main()
