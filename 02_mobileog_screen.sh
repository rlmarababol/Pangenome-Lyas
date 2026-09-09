#!/bin/bash
# 02_mobileog_screen.sh
# Mobile genetic element / horizontal-gene-transfer-associated gene
# screening (DIAMOND vs. mobileOG-db) -- Supplementary Methods S8.
# Requires predicted protein sequences (.faa) as input, e.g. from Bakta
# annotation output (00a_assembly_qc/06_annotate.sh).
# Usage: ./02_mobileog_screen.sh <proteins.faa> <mobileog_db.dmnd> <output.csv>

set -euo pipefail

PROTEINS="$1"
MOBILEOG_DB="$2"
OUTPUT="$3"

echo "=== Mobile genetic element screening (DIAMOND vs. mobileOG-db) ==="
diamond blastp --db "$MOBILEOG_DB" --query "$PROTEINS" \
  --outfmt 6 qseqid sseqid pident bitscore slen evalue qlen sstart send qstart qend \
  --evalue 1e-10 --out "$OUTPUT"

echo "Done. Result: $OUTPUT"
echo "NOTE: cross-reference hit IDs against the mobileOG-db metadata table to"
echo "recover gene name, major/minor category, and source database for each hit."
