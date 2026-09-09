#!/usr/bin/env python3
"""
Compute dN, dS, and dN/dS PER GENE FAMILY (Core and Shell), using the
same self-tested NG86 implementation already validated in
calculate_dnds_per_genome.py -- reused here verbatim, only the
aggregation level changes (per family, across all pairs of genomes
within that family's own PPanGGOLiN alignment, rather than per genome
vs. a single reference).

WHY THIS SCRIPT EXISTS: the original per-family attempt used
Biopython's experimental codonalign module, which failed on 100% of
pairs (147,906/147,906 errors) in this environment. That was fixed by
writing this self-contained NG86 implementation, but at that point the
analysis pivoted straight to the per-genome version and the per-family
file was never regenerated with the working code -- any
dnds_per_family.csv on disk from that period is invalid leftover data,
not a bug in code reading it.

Genome-wide (ratio of SUMMED counts, not average of per-pair ratios)
convention is used here too, same as the per-genome version, for the
same statistical stability reason.

Excludes: families flagged as likely misclustered paralogs, and
genome pairs within a family that are not codon-frame-safe (checked
per pair, not per whole family).
"""

import csv
import glob
import itertools
import math
import os
import re

from Bio.Data import CodonTable

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
SCREEN_CSV = os.path.join(BASE_DIR, "family_identity_screen.csv")
OUT_CSV = os.path.join(BASE_DIR, "dnds_per_family.csv")

GAP_CHARS = {"-", "."}
BASES = "ACGT"

_table = CodonTable.unambiguous_dna_by_id[11]
FORWARD_TABLE = dict(_table.forward_table)
STOP_CODONS = set(_table.stop_codons)


def translate_codon(codon):
    if codon in STOP_CODONS:
        return "*"
    return FORWARD_TABLE.get(codon)


def synonymous_nonsynonymous_sites(codon):
    aa0 = translate_codon(codon)
    if aa0 is None:
        return None
    S = 0.0
    for pos in range(3):
        syn_count = 0
        for b in BASES:
            if b == codon[pos]:
                continue
            new_codon = codon[:pos] + b + codon[pos + 1:]
            aa1 = translate_codon(new_codon)
            if aa1 is not None and aa1 == aa0:
                syn_count += 1
        S += syn_count / 3.0
    return S, 3.0 - S


def classify_step(codon_a, codon_b):
    aa_a, aa_b = translate_codon(codon_a), translate_codon(codon_b)
    if aa_a is None or aa_b is None:
        return None
    return 1 if aa_a == aa_b else 0


def pairwise_codon_diff(codon1, codon2):
    diff_positions = [i for i in range(3) if codon1[i] != codon2[i]]
    if not diff_positions:
        return 0.0, 0.0

    total_Sd, total_Nd, n_paths = 0.0, 0.0, 0
    for perm in itertools.permutations(diff_positions):
        current = list(codon1)
        path_Sd, path_Nd = 0.0, 0.0
        valid = True
        for pos in perm:
            prev_codon = "".join(current)
            current[pos] = codon2[pos]
            new_codon = "".join(current)
            syn = classify_step(prev_codon, new_codon)
            if syn is None:
                valid = False
                break
            if syn == 1:
                path_Sd += 1
            else:
                path_Nd += 1
        if valid:
            total_Sd += path_Sd
            total_Nd += path_Nd
            n_paths += 1

    if n_paths == 0:
        return None
    return total_Sd / n_paths, total_Nd / n_paths


def jc_correct(p):
    if p is None or p >= 0.75:
        return None
    return -0.75 * math.log(1 - (4.0 / 3.0) * p)


def accumulate_ng86(seq1, seq2):
    n_codons = len(seq1) // 3
    S_sites, N_sites, Sd, Nd = 0.0, 0.0, 0.0, 0.0

    for i in range(n_codons):
        c1 = seq1[i * 3:i * 3 + 3].upper()
        c2 = seq2[i * 3:i * 3 + 3].upper()

        if any(ch in GAP_CHARS for ch in c1 + c2):
            continue
        if any(b not in BASES for b in c1 + c2):
            continue
        if c1 in STOP_CODONS or c2 in STOP_CODONS:
            continue

        site1 = synonymous_nonsynonymous_sites(c1)
        site2 = synonymous_nonsynonymous_sites(c2)
        if site1 is None or site2 is None:
            continue

        S_sites += (site1[0] + site2[0]) / 2.0
        N_sites += (site1[1] + site2[1]) / 2.0

        if c1 != c2:
            diff = pairwise_codon_diff(c1, c2)
            if diff is None:
                continue
            Sd += diff[0]
            Nd += diff[1]

    return S_sites, N_sites, Sd, Nd


def self_test():
    S_sites, N_sites, Sd, Nd = accumulate_ng86("TTTAAA", "TTCAAA")
    print(f"Self-test 1 (synonymous-only pair): Sd={Sd} Nd={Nd}")
    assert Nd == 0 and Sd > 0, "self-test FAILED"
    S2, N2, Sd2, Nd2 = accumulate_ng86("AAA", "AAT")
    print(f"Self-test 2 (nonsynonymous pair): Sd={Sd2} Nd={Nd2}")
    assert Nd2 > 0, "self-test FAILED"
    print("Self-test PASSED.\n")


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


def is_frame_safe_pair(seq1, seq2):
    if len(seq1) == 0 or len(seq1) % 3 != 0 or len(seq1) != len(seq2):
        return False
    for s in (seq1, seq2):
        for run in re.findall(r"[-.]+", s):
            if len(run) % 3 != 0:
                return False
    return True


def load_flagged(screen_csv):
    flagged = set()
    with open(screen_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["flagged_misclustered"].strip().lower() == "true":
                flagged.add((row["category"], row["family_id"]))
    return flagged


def main():
    self_test()
    flagged = load_flagged(SCREEN_CSV)
    rows = []

    for category, folder in FOLDERS.items():
        files = sorted(glob.glob(os.path.join(folder, "*.aln")))
        print(f"{category}: processing {len(files)} families...")

        n_processed, n_skip_misclustered, n_skip_no_pairs = 0, 0, 0

        for idx, path in enumerate(files, 1):
            family_id = os.path.splitext(os.path.basename(path))[0]

            if (category, family_id) in flagged:
                n_skip_misclustered += 1
                continue

            seqs = parse_fasta_aln(path)
            if len(seqs) < 2:
                continue

            S_total, N_total, Sd_total, Nd_total, n_pairs = 0.0, 0.0, 0.0, 0.0, 0

            for s1, s2 in itertools.combinations(seqs, 2):
                if not is_frame_safe_pair(s1, s2):
                    continue
                S, N, Sd, Nd = accumulate_ng86(s1.upper(), s2.upper())
                S_total += S
                N_total += N
                Sd_total += Sd
                Nd_total += Nd
                n_pairs += 1

            if n_pairs == 0 or S_total == 0 or N_total == 0:
                n_skip_no_pairs += 1
                continue

            pS, pN = Sd_total / S_total, Nd_total / N_total
            dS, dN = jc_correct(pS), jc_correct(pN)
            if dS is None or dN is None:
                n_skip_no_pairs += 1
                continue
            dnds = (dN / dS) if dS > 0 else None

            rows.append({
                "category": category, "family_id": family_id,
                "n_seqs": len(seqs), "n_pairs_used": n_pairs,
                "mean_dN": round(dN, 6), "mean_dS": round(dS, 6),
                "dN_dS": round(dnds, 6) if dnds is not None else "",
            })
            n_processed += 1

            if idx % 500 == 0:
                print(f"  ...{category}: {idx}/{len(files)}")

        print(f"{category}: {n_processed} families processed | "
              f"{n_skip_misclustered} skipped (misclustered) | "
              f"{n_skip_no_pairs} skipped (no valid pairs)\n")

    with open(OUT_CSV, "w", newline="") as f:
        fieldnames = ["category", "family_id", "n_seqs", "n_pairs_used",
                      "mean_dN", "mean_dS", "dN_dS"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {OUT_CSV}")


if __name__ == "__main__":
    main()
