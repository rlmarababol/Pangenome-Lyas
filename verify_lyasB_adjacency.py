#!/usr/bin/env python3
"""
Verify physical adjacency of LyasB's genome-specific genes using real
GFF3 coordinates AND strand (not just locus-tag-number proximity),
then prepare data for a gene-cluster visualization.
"""

import csv
import os
import re

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
OUT_JSON = os.path.expanduser("~/Downloads/lyasB_specific_genes.json")

TARGET_GENES = [
    "ODLLIK_02095", "ODLLIK_09145", "ODLLIK_09690",
    "ODLLIK_12825", "ODLLIK_12830", "ODLLIK_12835", "ODLLIK_16155",
]


def main():
    path = os.path.join(GFF3_DIR, "LyasB.gff3")
    genes = []

    with open(path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = cols[8]
            id_match = re.search(r"ID=([^;]+)", attrs)
            product_match = re.search(r"product=([^;]+)", attrs)
            if not id_match or id_match.group(1) not in TARGET_GENES:
                continue
            genes.append({
                "gene_id": id_match.group(1),
                "contig": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "product": product_match.group(1) if product_match else "hypothetical protein",
            })

    genes.sort(key=lambda g: (g["contig"], g["start"]))

    print("=== All 7 LyasB genome-specific genes, real coordinates ===")
    for g in genes:
        print(f"  {g['gene_id']} | contig={g['contig']} | {g['start']}-{g['end']} | "
              f"strand={g['strand']} | {g['product']}")

    print("\n=== Adjacency check (consecutive genes on the SAME contig) ===")
    for i in range(len(genes) - 1):
        g1, g2 = genes[i], genes[i + 1]
        if g1["contig"] == g2["contig"]:
            gap = g2["start"] - g1["end"]
            same_strand = "SAME strand" if g1["strand"] == g2["strand"] else "DIFFERENT strand"
            print(f"  {g1['gene_id']} -> {g2['gene_id']}: gap = {gap}bp, {same_strand}")
        else:
            print(f"  {g1['gene_id']} -> {g2['gene_id']}: DIFFERENT CONTIGS ({g1['contig']} vs {g2['contig']})")

    import json
    with open(OUT_JSON, "w") as f:
        json.dump(genes, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")
    print("Upload this file to build the visualization.")


if __name__ == "__main__":
    main()
