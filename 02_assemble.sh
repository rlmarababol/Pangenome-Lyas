#!/bin/bash
# 02_assemble.sh
# De novo genome assembly (SPAdes, isolate mode) -- Supplementary Methods S8.
# Usage: ./02_assemble.sh <R1_paired.fastq.gz> <R2_paired.fastq.gz> <output_dir> <memory_gb>

set -euo pipefail

R1="$1"
R2="$2"
OUTDIR="$3"
MEM_GB="${4:-16}"

echo "=== De novo assembly (SPAdes, isolate mode) ==="
spades.py --isolate -1 "$R1" -2 "$R2" -o "$OUTDIR" -m "$MEM_GB"

echo "Done. Assembly: $OUTDIR/contigs.fasta"
