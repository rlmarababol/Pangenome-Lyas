#!/usr/bin/env python3
"""
Check whether existing per-family DNA alignments (from ppanggolin msa
--source dna) are codon-frame-safe: alignment length divisible by 3,
and every gap run a multiple of 3bp (so gaps never shift the reading
frame). This determines whether dN/dS can be computed directly from
these alignments, or whether a protein-guided codon realignment
(translate -> align protein -> back-translate/thread onto DNA, i.e.
the standard PAL2NAL approach) is needed first.
"""

import glob
import os
import re
import random

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
SAMPLE_SIZE = 200
GAP_CHARS = "-."


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


def check_family(path):
    seqs = parse_fasta_aln(path)
    if len(seqs) < 2:
        return None
    aln_len = len(seqs[0])
    if aln_len == 0 or any(len(s) != aln_len for s in seqs):
        return None

    length_ok = (aln_len % 3 == 0)

    # Check every gap run (contiguous run of gap chars in any sequence)
    # is a multiple of 3
    gap_runs_ok = True
    for s in seqs:
        for run in re.findall(rf"[{GAP_CHARS}]+", s):
            if len(run) % 3 != 0:
                gap_runs_ok = False
                break
        if not gap_runs_ok:
            break

    return {"length_ok": length_ok, "gap_runs_ok": gap_runs_ok}


def main():
    random.seed(42)
    for category, folder in FOLDERS.items():
        files = sorted(glob.glob(os.path.join(folder, "*.aln")))
        sample = random.sample(files, min(SAMPLE_SIZE, len(files)))

        n_checked = 0
        n_length_ok = 0
        n_gap_ok = 0
        n_both_ok = 0

        for path in sample:
            result = check_family(path)
            if result is None:
                continue
            n_checked += 1
            if result["length_ok"]:
                n_length_ok += 1
            if result["gap_runs_ok"]:
                n_gap_ok += 1
            if result["length_ok"] and result["gap_runs_ok"]:
                n_both_ok += 1

        print(f"=== {category} (sampled {n_checked} families) ===")
        print(f"  Alignment length divisible by 3: {n_length_ok}/{n_checked} "
              f"({100*n_length_ok/n_checked:.1f}%)")
        print(f"  All gap runs multiple of 3:       {n_gap_ok}/{n_checked} "
              f"({100*n_gap_ok/n_checked:.1f}%)")
        print(f"  Codon-frame-safe overall:          {n_both_ok}/{n_checked} "
              f"({100*n_both_ok/n_checked:.1f}%)")
        print()


if __name__ == "__main__":
    main()
