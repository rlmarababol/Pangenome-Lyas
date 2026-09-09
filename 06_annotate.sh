#!/bin/bash
# 06_annotate.sh
# Genome annotation (Bakta v1.12.1) -- Supplementary Methods S8.
# Usage: ./06_annotate.sh <assembly.fasta> <bakta_db_path> <output_dir>

set -euo pipefail

ASSEMBLY="$1"
BAKTA_DB="$2"
OUTDIR="$3"

echo "=== Genome annotation (Bakta v1.12.1) ==="
bakta --db "$BAKTA_DB" --output "$OUTDIR" "$ASSEMBLY"

echo "Done. Annotation outputs (GFF3, GenBank, JSON, etc.) in: $OUTDIR"
