#!/usr/bin/env python3
"""
Unified extraction for ALL 6 LyasA/LyasB genome-specific gene regions,
with ADAPTIVE window sizing:
  - If the contig is short (<=20kb total), show the WHOLE contig --
    no padding, "cut short" naturally to what's actually there.
  - If long, extend outward from the target gene(s) gene-by-gene in
    each direction, stopping at the first intergenic gap >3kb (a
    natural signal of leaving the locally relevant neighborhood), up
    to a hard cap of 20kb total width.

Covers: LyasA's EHHGDB_20800 (transposase), and LyasB's R-M cluster
(ODLLIK_12825/12830/12835 treated as one region) plus the 4 remaining
individual genes.
"""

import json
import os
import re

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
MAX_WINDOW = 20000
GAP_STOP_THRESHOLD = 3000
SHORT_CONTIG_THRESHOLD = 20000

# (region_name, genome, target_gene_ids) -- cluster genes grouped as one region
REGIONS = [
    ("LyasB_RM_cluster", "LyasB", ["ODLLIK_12825", "ODLLIK_12830", "ODLLIK_12835"]),
    ("LyasA_transposase", "LyasA", ["EHHGDB_20800"]),
    ("LyasB_02095", "LyasB", ["ODLLIK_02095"]),
    ("LyasB_09145", "LyasB", ["ODLLIK_09145"]),
    ("LyasB_09690", "LyasB", ["ODLLIK_09690"]),
    ("LyasB_16155", "LyasB", ["ODLLIK_16155"]),
]


def parse_gff3(genome):
    path = os.path.join(GFF3_DIR, f"{genome}.gff3")
    seq_lengths = {}
    all_cds = []
    with open(path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    seq_lengths[parts[1]] = int(parts[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = cols[8]
            id_match = re.search(r"ID=([^;]+)", attrs)
            product_match = re.search(r"product=([^;]+)", attrs)
            if not id_match:
                continue
            all_cds.append({
                "gene_id": id_match.group(1), "contig": cols[0],
                "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6],
                "product": product_match.group(1) if product_match else "hypothetical protein",
            })
    return seq_lengths, all_cds


def determine_window(target_genes, all_cds, contig_len):
    contig = target_genes[0]["contig"]
    contig_genes = sorted([g for g in all_cds if g["contig"] == contig], key=lambda g: g["start"])

    if contig_len <= SHORT_CONTIG_THRESHOLD:
        return 0, contig_len, "full contig (short enough)"

    target_ids = {g["gene_id"] for g in target_genes}
    target_min = min(g["start"] for g in target_genes)
    target_max = max(g["end"] for g in target_genes)

    # BUGFIX: previously checked each direction's extension against MAX_WINDOW
    # independently (relative to a fixed anchor), which allowed the TOTAL
    # combined span to reach up to ~2x MAX_WINDOW (confirmed: LyasB_16155
    # came out as a 40kb window instead of the intended 20kb cap). Now
    # checks the ACTUAL TOTAL width (right_bound - left_bound) at every
    # extension step, on whichever side is being extended.
    left_bound = target_min
    right_bound = target_max
    idx = next(i for i, g in enumerate(contig_genes) if g["gene_id"] in target_ids)

    left_idx = idx - 1
    right_idx = idx + 1
    while True:
        extended = False
        # try extending left
        if left_idx >= 0:
            candidate_gene = contig_genes[left_idx]
            gap = left_bound - candidate_gene["end"]
            candidate_left = candidate_gene["start"]
            if gap <= GAP_STOP_THRESHOLD and (right_bound - candidate_left) <= MAX_WINDOW:
                left_bound = candidate_left
                left_idx -= 1
                extended = True
        # try extending right
        if right_idx < len(contig_genes):
            candidate_gene = contig_genes[right_idx]
            gap = candidate_gene["start"] - right_bound
            candidate_right = candidate_gene["end"]
            if gap <= GAP_STOP_THRESHOLD and (candidate_right - left_bound) <= MAX_WINDOW:
                right_bound = candidate_right
                right_idx += 1
                extended = True
        if not extended:
            break

    window_start = max(0, left_bound - 200)
    window_end = min(contig_len, right_bound + 200)
    return window_start, window_end, f"adaptive (gap-stop, TOTAL capped at {MAX_WINDOW}bp)"


def main():
    results = {}
    genome_cache = {}

    for region_name, genome, target_ids in REGIONS:
        if genome not in genome_cache:
            genome_cache[genome] = parse_gff3(genome)
        seq_lengths, all_cds = genome_cache[genome]

        target_genes = [g for g in all_cds if g["gene_id"] in target_ids]
        if not target_genes:
            print(f"WARNING: none of {target_ids} found in {genome}\n")
            continue

        contig = target_genes[0]["contig"]
        contig_len = seq_lengths.get(contig)
        if contig_len is None:
            print(f"WARNING: no length for {contig}\n")
            continue

        window_start, window_end, sizing = determine_window(target_genes, all_cds, contig_len)

        genes_in_window = [
            g for g in all_cds
            if g["contig"] == contig and g["end"] >= window_start and g["start"] <= window_end
        ]
        genes_in_window.sort(key=lambda g: g["start"])

        print(f"=== {region_name} ({genome}, {contig}) ===")
        print(f"  Contig length: {contig_len}bp | sizing: {sizing}")
        print(f"  Window: {window_start}-{window_end} ({window_end - window_start}bp)")
        print(f"  Genes shown: {len(genes_in_window)}\n")

        results[region_name] = {
            "genome": genome, "contig": contig, "contig_length": contig_len,
            "window_start": window_start, "window_end": window_end,
            "target_genes": target_ids,
            "genes": genes_in_window,
        }

    out_json = os.path.expanduser("~/Downloads/all_gene_contexts.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_json}")
    print("Upload this file to build all 6 visualizations consistently.")


if __name__ == "__main__":
    main()
