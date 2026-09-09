#!/bin/bash
# 03_assembly_qc.sh
# Assembly contiguity and completeness assessment (QUAST + CheckM2) --
# Supplementary Methods S8; completeness values reported in Table 1
# were generated using CheckM2 v1.1.0.
# Usage: ./03_assembly_qc.sh <contigs.fasta> <output_dir> <threads>

set -euo pipefail

CONTIGS="$1"
OUTDIR="$2"
THREADS="${3:-4}"

mkdir -p "$OUTDIR"

echo "=== Assembly contiguity (QUAST) ==="
quast.py "$CONTIGS" -o "$OUTDIR/quast_out"

echo "=== Genome completeness (CheckM2 v1.1.0) ==="
checkm2 predict --threads "$THREADS" --input "$CONTIGS" --output-directory "$OUTDIR/checkm2_out"

echo "Done."
echo "  Contiguity report: $OUTDIR/quast_out/report.txt"
echo "  Completeness report: $OUTDIR/checkm2_out/quality_report.tsv"
