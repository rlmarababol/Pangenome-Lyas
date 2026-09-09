#!/usr/bin/env python3
"""
Classify large-effect mutations per genome, relative to reference
genome LyasI, using the same PPanGGOLiN Core/Shell family alignments
as the per-genome SNP/indel pipeline (count_snp_indel_per_genome.py).

Categories:
  - Frameshift indel: an indel (contiguous gap run) between LyasI and
    a query genome whose length is NOT a multiple of 3bp.
  - In-frame indel: same, but length IS a multiple of 3bp (whole-codon
    insertion/deletion -- reading frame preserved).
  - Stop-gain: a codon where LyasI has a sense codon and the query has
    a stop codon (premature termination in the query, relative to ref).
  - Stop-loss: the reverse -- LyasI has a stop codon, query does not
    (read-through/extension in the query).

SCOPE (stated explicitly, not a silent limitation):
  Stop-gain/stop-loss are only classified for gene-pairs (one genome vs
  LyasI, within one family) that have ZERO indels between them. This
  keeps both sequences the same length and correctly in-frame, so
  codon-by-codon comparison is simple and fully correct. After a real
  frameshift, everything downstream translates in a different frame;
  classifying stop codons there would require full frame-tracking,
  which this script does not attempt. Gene-pairs with any indel are
  excluded from stop-codon classification (the indel itself is still
  reported as frameshift/in-frame).

Excludes: families flagged as likely misclustered paralogs, and
families where LyasI is not present.
"""

import csv
import glob
import os

from Bio.Data import CodonTable

BASE_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results")
FOLDERS = {
    "Core":  os.path.join(BASE_DIR, "msa_core", "msa_core_dna"),
    "Shell": os.path.join(BASE_DIR, "msa_shell", "msa_shell_dna"),
}
SCREEN_CSV = os.path.join(BASE_DIR, "family_identity_screen.csv")
SUMMARY_CSV = os.path.join(BASE_DIR, "large_effect_summary.csv")
DETAIL_CSV = os.path.join(BASE_DIR, "large_effect_detail.csv")

REFERENCE = "LyasI"
GAP_CHARS = {"-", "."}
BASES = "ACGT"

_table = CodonTable.unambiguous_dna_by_id[11]
FORWARD_TABLE = dict(_table.forward_table)
STOP_CODONS = set(_table.stop_codons)


def translate_codon(codon):
    if codon in STOP_CODONS:
        return "*"
    return FORWARD_TABLE.get(codon)


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


def load_flagged(screen_csv):
    flagged = set()
    with open(screen_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["flagged_misclustered"].strip().lower() == "true":
                flagged.add((row["category"], row["family_id"]))
    return flagged


def classify_pair(ref_seq, query_seq):
    """
    ref_seq, query_seq: raw alignment rows (with gaps) for LyasI and
    one query genome, same family alignment, same length.
    Returns dict with 'indels' (list of {length, ref_pos, type}) and,
    if the pair has zero indels, 'stop_gain' / 'stop_loss' (lists of
    codon positions).
    """
    ref_seq, query_seq = ref_seq.upper(), query_seq.upper()
    cleaned = [(r, q) for r, q in zip(ref_seq, query_seq)
               if not (r in GAP_CHARS and q in GAP_CHARS)]

    indels = []
    ref_pos = 0
    run_len = 0
    run_start = None
    run_type = None

    for r, q in cleaned:
        r_gap, q_gap = r in GAP_CHARS, q in GAP_CHARS
        if r_gap != q_gap:
            if run_len == 0:
                run_start = ref_pos
                run_type = "insertion_in_query" if r_gap else "deletion_in_query"
            run_len += 1
        else:
            if run_len > 0:
                indels.append({"length": run_len, "ref_pos": run_start, "type": run_type})
                run_len = 0
        if not r_gap:
            ref_pos += 1
    if run_len > 0:
        indels.append({"length": run_len, "ref_pos": run_start, "type": run_type})

    result = {"indels": indels, "stop_gain": [], "stop_loss": []}

    if len(indels) == 0:
        ref_ungapped = "".join(r for r, q in cleaned)
        query_ungapped = "".join(q for r, q in cleaned)
        n_codons = len(ref_ungapped) // 3
        for c in range(n_codons):
            rc = ref_ungapped[c * 3:c * 3 + 3]
            qc = query_ungapped[c * 3:c * 3 + 3]
            if rc == qc:
                continue
            if any(b not in BASES for b in rc + qc):
                continue
            ref_is_stop = rc in STOP_CODONS
            query_is_stop = qc in STOP_CODONS
            if not ref_is_stop and query_is_stop:
                result["stop_gain"].append(c)
            elif ref_is_stop and not query_is_stop:
                result["stop_loss"].append(c)

    return result


def main():
    flagged = load_flagged(SCREEN_CSV)

    # (genome, category) -> counts
    summary = {}
    detail_rows = []

    for category, folder in FOLDERS.items():
        files = sorted(glob.glob(os.path.join(folder, "*.aln")))
        print(f"{category}: scanning {len(files)} families...")

        n_used = 0
        n_no_reference = 0
        n_misclustered = 0

        for idx, path in enumerate(files, 1):
            family_id = os.path.splitext(os.path.basename(path))[0]

            if (category, family_id) in flagged:
                n_misclustered += 1
                continue

            seqs = parse_fasta_aln(path)
            if REFERENCE not in seqs:
                n_no_reference += 1
                continue

            ref_seq = seqs[REFERENCE]
            n_used += 1

            for genome, seq in seqs.items():
                if genome == REFERENCE:
                    continue

                result = classify_pair(ref_seq, seq)
                key = (genome, category)
                if key not in summary:
                    summary[key] = {
                        "n_families_compared": 0, "n_frameshift": 0, "n_inframe_indel": 0,
                        "n_stop_gain": 0, "n_stop_loss": 0,
                        "n_families_skipped_for_stopcodon": 0,
                    }
                summary[key]["n_families_compared"] += 1

                for indel in result["indels"]:
                    is_frameshift = (indel["length"] % 3 != 0)
                    event_type = "frameshift" if is_frameshift else "in_frame_indel"
                    summary[key]["n_frameshift" if is_frameshift else "n_inframe_indel"] += 1
                    detail_rows.append({
                        "genome": genome, "category": category, "family_id": family_id,
                        "event_type": event_type, "position": indel["ref_pos"],
                        "length_bp": indel["length"],
                    })

                if len(result["indels"]) > 0:
                    summary[key]["n_families_skipped_for_stopcodon"] += 1

                for codon_pos in result["stop_gain"]:
                    summary[key]["n_stop_gain"] += 1
                    detail_rows.append({
                        "genome": genome, "category": category, "family_id": family_id,
                        "event_type": "stop_gain", "position": codon_pos, "length_bp": "",
                    })
                for codon_pos in result["stop_loss"]:
                    summary[key]["n_stop_loss"] += 1
                    detail_rows.append({
                        "genome": genome, "category": category, "family_id": family_id,
                        "event_type": "stop_loss", "position": codon_pos, "length_bp": "",
                    })

            if idx % 500 == 0:
                print(f"  ...{category}: {idx}/{len(files)}")

        print(f"{category}: {n_used} families used | "
              f"{n_no_reference} skipped (no {REFERENCE}) | "
              f"{n_misclustered} skipped (misclustered)\n")

    summary_rows = []
    for (genome, category), c in summary.items():
        summary_rows.append({
            "genome": genome, "category": category,
            "n_families_compared": c["n_families_compared"],
            "n_families_skipped_for_stopcodon": c["n_families_skipped_for_stopcodon"],
            "n_stop_gain": c["n_stop_gain"],
            "n_stop_loss": c["n_stop_loss"],
            "n_frameshift": c["n_frameshift"],
            "n_inframe_indel": c["n_inframe_indel"],
        })
    summary_rows.sort(key=lambda r: (r["category"], r["genome"]))

    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    with open(DETAIL_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["genome", "category", "family_id",
                                                 "event_type", "position", "length_bp"])
        writer.writeheader()
        writer.writerows(detail_rows)

    print("=== Per-genome summary ===")
    for row in summary_rows:
        print(f"  {row['genome']:8s} {row['category']:6s} | "
              f"stop_gain={row['n_stop_gain']:3d}  stop_loss={row['n_stop_loss']:3d}  "
              f"frameshift={row['n_frameshift']:3d}  in_frame_indel={row['n_inframe_indel']:3d}")

    print(f"\nSaved: {SUMMARY_CSV}")
    print(f"Saved: {DETAIL_CSV}")


if __name__ == "__main__":
    main()
