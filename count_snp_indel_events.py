#!/usr/bin/env python3
"""
Count SNP and indel events per gene family alignment (Core vs Shell),
normalized to events per kb of alignment length. One row per gene
family = one data point for the downstream R plot.

v2 change (methodological fix): with up to 10 draft genome assemblies
compared simultaneously, a single genome's sequencing/assembly error at
any column is enough to flag it as a "SNP" or "indel" under a naive
"do all sequences agree" rule -- and with more genomes being compared
at once, the odds that at least one has a stray error at any given
position rise accordingly. To guard against this:
  - Ambiguous IUPAC codes (N, R, Y, S, W, K, M, B, D, H, V) are treated
    as missing data, never as evidence of a distinct allele.
  - The minority state (the rarer base, or the rarer gap/non-gap split)
    must be supported by at least MIN_MINOR_COUNT genomes -- not just
    one -- to be called a real SNP or indel.
This is a genuine sensitivity/specificity trade-off (it will miss true
singleton variants carried by only one genome) in exchange for filtering
out assembly noise -- the standard, defensible choice for unpolished
draft assemblies. Adjust MIN_MINOR_COUNT below to change that trade-off.
"""

import csv
import glob
import os

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
OUT_CSV = os.path.join(BASE_DIR, "snp_indel_events.csv")

GAP_CHARS = {"-", "."}
VALID_BASES = {"A", "C", "G", "T"}   # ambiguous/N codes excluded, treated as missing
MIN_MINOR_FRACTION = 0.20             # minor state must represent >=20% of genomes
                                       # present, scales correctly regardless of family
                                       # size (fixes the old rule's blind spot where any
                                       # 2-genome family could never register a variant)


def parse_fasta_aln(path):
    """Read a simple multi-FASTA alignment file into a list of sequences."""
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


def count_events(seqs, min_minor_fraction=MIN_MINOR_FRACTION):
    n = len(seqs)
    if n < 2:
        return None

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        return None
    aln_len = lengths.pop()
    if aln_len == 0:
        return None

    snp_count = 0
    indel_col = [False] * aln_len

    for i in range(aln_len):
        col = [s[i].upper() for s in seqs]
        gap_count = sum(1 for c in col if c in GAP_CHARS)
        nongap_count = n - gap_count

        if gap_count > 0 and nongap_count > 0:
            # Candidate indel column: minority side (gap or non-gap,
            # whichever is rarer) must clear the proportional threshold
            minor_count = min(gap_count, nongap_count)
            if minor_count / n >= min_minor_fraction:
                indel_col[i] = True
        elif gap_count == 0:
            # Candidate SNP column: only confident A/C/G/T calls count;
            # N/ambiguous codes are treated as missing, not as an allele.
            confident = [c for c in col if c in VALID_BASES]
            n_confident = len(confident)
            if n_confident < 2:
                continue  # not enough confident calls to evaluate this site
            counts = {}
            for b in confident:
                counts[b] = counts.get(b, 0) + 1
            if len(counts) > 1:
                minor_count = sorted(counts.values())[-2]
                if minor_count / n_confident >= min_minor_fraction:
                    snp_count += 1

    indel_events = 0
    prev = False
    for v in indel_col:
        if v and not prev:
            indel_events += 1
        prev = v

    total_events = snp_count + indel_events
    events_per_kb = total_events / (aln_len / 1000)

    return {
        "n_seqs": n,
        "aln_len": aln_len,
        "snp": snp_count,
        "indel_events": indel_events,
        "total_events": total_events,
        "events_per_kb": events_per_kb,
    }


def main():
    rows = []
    skipped = {"Core": 0, "Shell": 0}
    processed = {"Core": 0, "Shell": 0}

    for category, folder in FOLDERS.items():
        files = sorted(glob.glob(os.path.join(folder, "*.aln")))
        print(f"{category}: found {len(files)} alignment files in {folder}")

        for idx, path in enumerate(files, 1):
            family_id = os.path.splitext(os.path.basename(path))[0]
            seqs = parse_fasta_aln(path)
            result = count_events(seqs)

            if result is None:
                skipped[category] += 1
                continue

            rows.append({
                "category": category,
                "family_id": family_id,
                **result,
            })
            processed[category] += 1

            if idx % 500 == 0:
                print(f"  ...{category}: {idx}/{len(files)} files processed")

    with open(OUT_CSV, "w", newline="") as f:
        fieldnames = ["category", "family_id", "n_seqs", "aln_len",
                      "snp", "indel_events", "total_events", "events_per_kb"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("\n=== Summary ===")
    for cat in FOLDERS:
        print(f"{cat}: {processed[cat]} families used, {skipped[cat]} skipped "
              f"(fewer than 2 sequences present, or malformed alignment)")
    print(f"\nMIN_MINOR_FRACTION = {MIN_MINOR_FRACTION} (variant must represent "
          f"this fraction of genomes present to be counted)")
    print(f"Saved: {OUT_CSV}")


if __name__ == "__main__":
    main()
