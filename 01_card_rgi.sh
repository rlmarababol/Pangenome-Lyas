#!/bin/bash
# 01_card_rgi.sh
# Antibiotic resistance gene screening (RGI against CARD, strict
# cut-off) -- Supplementary Methods S8.
# Usage: ./01_card_rgi.sh <contigs.fasta> <output_dir>

set -euo pipefail

CONTIGS="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"

echo "=== Antibiotic resistance gene screening (RGI / CARD, strict cut-off) ==="
rgi main --input_sequence "$CONTIGS" --output_file "$OUTDIR/rgi_out" --local -a DIAMOND --clean

echo "Done. Result: ${OUTDIR}/rgi_out.txt"
echo "NOTE: interpret low-identity hits (<~50%) with caution -- these indicate"
echo "distant sequence relatedness, not confirmed functional resistance."
