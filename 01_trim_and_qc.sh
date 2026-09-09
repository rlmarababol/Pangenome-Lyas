#!/bin/bash
# 01_trim_and_qc.sh
# Read trimming and quality control -- Supplementary Methods S8.
# Usage: ./01_trim_and_qc.sh <sample_prefix> <R1.fastq.gz> <R2.fastq.gz> <output_dir> <threads> <heap_gb>

set -euo pipefail

SAMPLE="$1"
R1="$2"
R2="$3"
OUTDIR="$4"
THREADS="${5:-4}"
HEAP_GB="${6:-4}"

mkdir -p "$OUTDIR"

echo "=== Trimming (Trimmomatic, min Phred Q30) ==="
trimmomatic PE -threads "$THREADS" -Xmx"${HEAP_GB}"g \
  "$R1" "$R2" \
  "$OUTDIR/${SAMPLE}_R1_paired.fastq.gz" "$OUTDIR/${SAMPLE}_R1_unpaired.fastq.gz" \
  "$OUTDIR/${SAMPLE}_R2_paired.fastq.gz" "$OUTDIR/${SAMPLE}_R2_unpaired.fastq.gz" \
  LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:50

echo "=== Read QC (FastQC) ==="
fastqc "$OUTDIR/${SAMPLE}_R1_paired.fastq.gz" "$OUTDIR/${SAMPLE}_R2_paired.fastq.gz" -o "$OUTDIR"

echo "Done. Trimmed, paired reads: $OUTDIR/${SAMPLE}_R1_paired.fastq.gz / ${SAMPLE}_R2_paired.fastq.gz"
