#!/usr/bin/env python3
"""
Identify genome-specific (prevalence=1) genes for LyasA and LyasB
specifically, with product annotations.
"""

import csv
import os
import re

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
SYN_DIR = os.path.join(BASE_DIR, "synteny")
GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
RTAB = os.path.join(BASE_DIR, "matrix_export", "gene_presence_absence.Rtab")
OUT_CSV = os.path.join(BASE_DIR, "genome_specific_lyasA_lyasB.csv")

TARGET_GENOMES = {"LyasA", "LyasB"}


def load_category_map():
    core_ids = set(open(os.path.join(BASE_DIR, "matrix_export/partitions/persistent.txt")).read().split())
    shell_ids = set(open(os.path.join(BASE_DIR, "matrix_export/partitions/shell.txt")).read().split())
    cloud_ids = set(open(os.path.join(BASE_DIR, "matrix_export/partitions/cloud.txt")).read().split())
    return core_ids, shell_ids, cloud_ids


def load_gene_to_family():
    """gene_id -> correct (.aln-matching) family_id, from genes.tsv."""
    mapping = {}
    with open(os.path.join(SYN_DIR, "genes.tsv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["family_id"]:
                mapping[row["gene_id"]] = row["family_id"]
    return mapping


def find_specific_families():
    """Families present in exactly 1 genome, restricted to LyasA/LyasB."""
    results = []
    with open(RTAB) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        genome_cols = header[1:]
        for row in reader:
            rtab_family_id = row[0]
            presence = row[1:]
            present_genomes = [g for g, p in zip(genome_cols, presence) if p == "1"]
            if len(present_genomes) == 1 and present_genomes[0] in TARGET_GENOMES:
                results.append((rtab_family_id, present_genomes[0]))
    return results


def extract_products_and_genes():
    """genome -> gene_id -> (product, start, end) parsed fresh from GFF3."""
    gene_info = {}
    for genome_file in os.listdir(GFF3_DIR):
        if not genome_file.endswith(".gff3"):
            continue
        genome_name = genome_file.replace(".gff3", "")
        if genome_name not in TARGET_GENOMES:
            continue
        gene_info[genome_name] = {}
        with open(os.path.join(GFF3_DIR, genome_file)) as f:
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
                if id_match:
                    product = product_match.group(1) if product_match else "hypothetical protein"
                    gene_info[genome_name][id_match.group(1)] = {
                        "product": product, "contig": cols[0],
                        "start": cols[3], "end": cols[4],
                    }
    return gene_info


def main():
    core_ids, shell_ids, cloud_ids = load_category_map()
    gene_to_family = load_gene_to_family()
    gene_info = extract_products_and_genes()

    # reverse: family_id -> the specific gene_id in that genome
    family_to_gene = {}
    for gene_id, family_id in gene_to_family.items():
        family_to_gene.setdefault(family_id, []).append(gene_id)

    specific_families = find_specific_families()
    print(f"Genome-specific families found for LyasA/LyasB: {len(specific_families)}")

    rows = []
    for rtab_family_id, genome in specific_families:
        family_id = gene_to_family.get(rtab_family_id, rtab_family_id)

        if family_id in core_ids:
            category = "Core"
        elif family_id in shell_ids:
            category = "Shell"
        elif family_id in cloud_ids:
            category = "Cloud"
        else:
            category = "Unknown"

        # find the actual gene_id belonging to THIS genome within the family
        candidates = family_to_gene.get(family_id, [family_id])
        this_genome_gene = None
        for g in candidates:
            if g in gene_info.get(genome, {}):
                this_genome_gene = g
                break
        if this_genome_gene is None and rtab_family_id in gene_info.get(genome, {}):
            this_genome_gene = rtab_family_id

        info = gene_info.get(genome, {}).get(this_genome_gene, {})
        rows.append({
            "genome": genome, "family_id": family_id, "gene_id": this_genome_gene or "?",
            "category": category, "product": info.get("product", "?"),
            "contig": info.get("contig", "?"), "start": info.get("start", "?"), "end": info.get("end", "?"),
        })

    rows.sort(key=lambda r: (r["genome"], r["gene_id"] or ""))

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {OUT_CSV}\n")

    for genome in sorted(TARGET_GENOMES):
        genome_rows = [r for r in rows if r["genome"] == genome]
        print(f"=== {genome}: {len(genome_rows)} genome-specific genes ===")
        for r in genome_rows:
            print(f"  {r['gene_id']} ({r['category']}) -- {r['product']}")
        print()


if __name__ == "__main__":
    main()
