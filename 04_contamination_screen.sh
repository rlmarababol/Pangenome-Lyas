#!/bin/bash
# 04_contamination_screen.sh
# Taxonomic contamination screening of assembled contigs (BLASTN vs. nt)
# -- Supplementary Methods S8. Contigs with top hits to non-target taxa
# should be flagged for manual review and removal prior to downstream
# analysis.
# Usage: ./04_contamination_screen.sh <contigs.fasta> <blast_nt_db_path> <output_dir>

set -euo pipefail

CONTIGS="$1"
NT_DB="$2"
OUTDIR="$3"

mkdir -p "$OUTDIR"

echo "=== Contamination screening (BLASTN vs. nt) ==="
blastn -query "$CONTIGS" -db "$NT_DB" -outfmt 6 -max_target_seqs 5 -evalue 1e-10 \
  -out "$OUTDIR/contamination_screen.tsv"

echo "Done. Review top hits per contig in: $OUTDIR/contamination_screen.tsv"
echo "Contigs whose top hit is not the target genus/species should be flagged for removal."
