#!/usr/bin/env python3
"""
Extract gene_id -> product annotation from each genome's GFF3 (CDS
features, "product=" attribute), so the top gene-family bar chart can
show actual gene names rather than opaque locus tags.
"""

import csv
import os
import re

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
OUT_PATH = os.path.expanduser(
    "~/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/gene_products.csv")

GENOMES = ["LyasA", "LyasB", "LyasC", "LyasD", "LyasE",
           "LyasF", "LyasG", "LyasH", "LyasI", "LyasJ"]


def parse_products(path, genome_name):
    rows = []
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
            if not id_match:
                continue
            gene_id = id_match.group(1)
            product = product_match.group(1) if product_match else "hypothetical protein"
            rows.append({"genome": genome_name, "gene_id": gene_id, "product": product})
    return rows


def main():
    all_rows = []
    for genome in GENOMES:
        path = os.path.join(GFF3_DIR, f"{genome}.gff3")
        rows = parse_products(path, genome)
        all_rows.extend(rows)
        print(f"{genome}: {len(rows)} gene products extracted")

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["genome", "gene_id", "product"])
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nTotal: {len(all_rows)} gene product annotations")
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
