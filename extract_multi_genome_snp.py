#!/usr/bin/env python3
"""
Extract a multi-genome alignment window centered on the single SNP
event in a chosen family, for the same style of multi-genome
alignment figure built for the indel example.

Run locally, then upload the resulting snp_alignment_data.json.
"""

import json
import os

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
GAP_CHARS = {"-", "."}

FAMILY_ID = "EHHGDB_02965"   # the C/E-matching example; change to visualize a different family
CATEGORY = "Core"

TREE_ORDER = ["LyasH", "LyasD", "LyasG", "LyasE", "LyasC",
              "LyasJ", "LyasF", "LyasB", "LyasA", "LyasI"]


def translate_codon(c):
    return CODON_TABLE.get(c.upper(), "X")


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


def main():
    base_dir = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
    msa_dir = "msa_core" if CATEGORY == "Core" else "msa_shell"
    aln_path = os.path.join(base_dir, msa_dir, f"{msa_dir}_dna", f"{FAMILY_ID}.aln")

    seqs = parse_fasta_aln(aln_path)
    aln_len = len(next(iter(seqs.values())))

    # Re-find the SNP column directly (same logic as the finder script)
    names = list(seqs.keys())
    snp_col = None
    minority_genomes = []
    for col in range(aln_len):
        col_chars = {n: seqs[n][col].upper() for n in names}
        if any(c in GAP_CHARS for c in col_chars.values()):
            continue
        distinct = set(col_chars.values())
        if len(distinct) > 1:
            snp_col = col
            base_to_genomes = {}
            for n, b in col_chars.items():
                base_to_genomes.setdefault(b, []).append(n)
            sorted_alleles = sorted(base_to_genomes.items(), key=lambda x: len(x[1]))
            minority_genomes = sorted_alleles[0][1]
            break

    print(f"SNP event: column {snp_col}, minority allele in: {minority_genomes}")

    # Window: codon-aligned, centered on the SNP, ~90bp total
    window_start = max(0, snp_col - 45)
    window_start -= window_start % 3
    window_end = min(aln_len, window_start + 90)
    print(f"Window: {window_start}-{window_end} ({window_end - window_start}bp)")

    rows = []
    for genome in TREE_ORDER:
        if genome not in seqs:
            continue
        seq = seqs[genome][window_start:window_end].upper()
        codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        aa = []
        for c in codons:
            if any(ch in GAP_CHARS for ch in c):
                aa.append("-")
            else:
                aa.append(translate_codon(c))
        rows.append({"genome": genome, "codons": codons, "aa": aa,
                     "has_minority_snp": genome in minority_genomes})

    data = {
        "family_id": FAMILY_ID, "category": CATEGORY,
        "snp_column": snp_col,
        "window_start": window_start, "window_end": window_end,
        "minority_genomes": minority_genomes,
        "rows": rows,
    }

    out_path = os.path.expanduser("~/Downloads/snp_alignment_data.json")
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nSaved: {out_path}")
    print("Upload this file back to continue building the visualization.")


if __name__ == "__main__":
    main()
