#!/usr/bin/env python3
"""
Extract and translate FULL protein sequences for specified genes,
looked up directly by gene ID from the correct genome's GFF3 (no
manual coordinate entry needed) -- ready to paste into NCBI BLASTP or
InterProScan.
"""

import re

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

GFF3_DIR = "/Users/ramjuliusmarababol/Downloads/GFF3"

# (gene_id, genome) -- all LyasA hypothetical genome-specific genes,
# plus LyasB's second hypothetical (12835 already done separately, but
# included again here for a complete, single combined output)
TARGETS = [
    ("EHHGDB_07575", "LyasA"),
    ("EHHGDB_12085", "LyasA"),
    ("EHHGDB_13870", "LyasA"),
    ("EHHGDB_18250", "LyasA"),
    ("EHHGDB_19375", "LyasA"),
    ("EHHGDB_19380", "LyasA"),
    ("EHHGDB_20790", "LyasA"),
    ("ODLLIK_12835", "LyasB"),
    ("ODLLIK_16155", "LyasB"),
]


def revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def translate(nt_seq):
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i+3].upper()
        aa.append(CODON_TABLE.get(codon, "X"))
    if aa and aa[-1] == "*":
        aa = aa[:-1]
    return "".join(aa)


def find_gene_coords(gff3_path, gene_id):
    """Search the GFF3 directly for a CDS matching this gene ID."""
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            if f"ID={gene_id}" in cols[8] or f"ID={gene_id};" in cols[8]:
                return {"contig": cols[0], "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6]}
    return None


def load_contig_sequences(gff3_path):
    sequences = {}
    in_fasta = False
    current_id = None
    current_seq = []
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if not in_fasta:
                continue
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences


def main():
    genome_cache = {}
    results = []

    for gene_id, genome in TARGETS:
        if genome not in genome_cache:
            path = f"{GFF3_DIR}/{genome}.gff3"
            genome_cache[genome] = {
                "coords_source": path,
                "contigs": load_contig_sequences(path),
            }

        coords = find_gene_coords(genome_cache[genome]["coords_source"], gene_id)
        if coords is None:
            print(f"WARNING: {gene_id} not found in {genome}.gff3\n")
            continue

        contig_seq = genome_cache[genome]["contigs"].get(coords["contig"])
        if not contig_seq:
            print(f"WARNING: contig {coords['contig']} not found for {gene_id}\n")
            continue

        nt = contig_seq[coords["start"] - 1:coords["end"]]
        if coords["strand"] == "-":
            nt = revcomp(nt)
        aa = translate(nt)

        header = f">{gene_id} [{genome}] {coords['contig']}:{coords['start']}-{coords['end']}({coords['strand']}) length={len(aa)}aa"
        print(header)
        for i in range(0, len(aa), 60):
            print(aa[i:i+60])
        print()

        results.append(header)

    print(f"\n=== {len(results)} / {len(TARGETS)} sequences extracted successfully ===")


if __name__ == "__main__":
    main()
