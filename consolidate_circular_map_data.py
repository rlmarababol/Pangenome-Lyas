#!/usr/bin/env python3
"""
Consolidate contig lengths + all annotation hits (CARD, mobileOG,
CRISPR) into one clean JSON per genome, for a Proksee/CGView-style
circular genome map (contigs concatenated sequentially, since these
are draft assemblies, not closed circular chromosomes).
"""

import csv
import json
import os
import re

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"
UPLOADS_DIR = os.path.expanduser("~/Downloads/Other common")

FILES = {
    "LyasA": {
        "card": "LyasA - CARD output.txt",
        "mobileog": "LyasA - cgview.fa.mobileOG.Alignment.Out.csv",
        "crispr": "LyasA - CRISPR-Cas_summary.tsv",
    },
    "LyasB": {
        "card": "LyasB - CARD output.txt",
        "mobileog": "LyasB - cgview.fa.mobileOG.Alignment.Out.csv",
        "crispr": None,  # not available
    },
}


def get_contig_lengths(genome):
    """From GFF3 ##sequence-region declarations. BUGFIX: previously
    broke on the first non-'#' line, incorrectly assuming all
    ##sequence-region declarations are clustered at the file's top --
    Bakta actually interspenses them per-contig alongside that
    contig's own features. Now reads the whole file (stopping only at
    ##FASTA, matching the established, validated pattern used
    elsewhere in this project) so all contigs are found regardless of
    where their declaration sits relative to feature lines."""
    path = os.path.join(GFF3_DIR, f"{genome}.gff3")
    lengths = {}
    with open(path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    lengths[parts[1]] = int(parts[3])
    return lengths


def node_name_to_contig(node_name, contig_lengths):
    """CARD/mobileOG use SPAdes NODE_ names; GFF3 uses contig_N. Match
    by LENGTH (embedded in the NODE name) since that's the only
    reliable shared identifier between the two naming systems."""
    m = re.search(r"length_(\d+)", node_name)
    if not m:
        return None
    length = int(m.group(1))
    for contig, clen in contig_lengths.items():
        if clen == length:
            return contig
    return None


def parse_card(path):
    hits = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hits.append({
                "node": row["Contig"], "start": int(row["Start"]), "end": int(row["Stop"]),
                "gene": row["Best_Hit_ARO"], "identity": float(row["Best_Identities"]),
                "drug_class": row["Drug Class"],
            })
    return hits


def parse_mobileog(path):
    hits = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            node = row["Contig/ORF Name"]
            hits.append({
                "node": node, "start": int(row["ORF_Start"]), "end": int(row["ORF_End"]),
                "gene": row["Gene Name"], "identity": float(row["Pident"]),
                "category": row["Major mobileOG Category"],
            })
    return hits


def parse_crispr(path):
    hits = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Nb CRISPRs"] and int(row["Nb CRISPRs"]) > 0:
                m = re.search(r"\[(\d+);(\d+)\]", row["CRISPR array(s)"])
                if m:
                    hits.append({"node": row["Sequence(s)"], "start": int(m.group(1)), "end": int(m.group(2))})
    return hits


def main():
    results = {}
    for genome, files in FILES.items():
        contig_lengths = get_contig_lengths(genome)

        card_hits = parse_card(os.path.join(UPLOADS_DIR, files["card"]))
        mobileog_hits = parse_mobileog(os.path.join(UPLOADS_DIR, files["mobileog"]))
        crispr_hits = parse_crispr(os.path.join(UPLOADS_DIR, files["crispr"])) if files["crispr"] else []

        def attach_contig(hits):
            for h in hits:
                h["contig"] = node_name_to_contig(h["node"], contig_lengths)
            return [h for h in hits if h["contig"] is not None]

        card_hits = attach_contig(card_hits)
        mobileog_hits = attach_contig(mobileog_hits)
        crispr_hits = attach_contig(crispr_hits)

        results[genome] = {
            "contig_lengths": contig_lengths,
            "card_hits": card_hits,
            "mobileog_hits": mobileog_hits,
            "crispr_hits": crispr_hits,
        }

        print(f"{genome}: {len(contig_lengths)} contigs, {len(card_hits)} CARD hits, "
              f"{len(mobileog_hits)} mobileOG hits matched to contigs, {len(crispr_hits)} CRISPR hits")

    out_path = os.path.expanduser("~/Downloads/circular_map_data.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")
    print("Upload this file to build the circular genome map.")


if __name__ == "__main__":
    main()
