#!/bin/bash
set -euo pipefail

# Usage: ./run_coverage_mapping.sh <sample_name> <r1_fastq> <r2_fastq> <assembly_fasta>
# Example (1c/LyasA):
#   ./run_coverage_mapping.sh LyasA \
#     ~/Downloads/macrogen_1c/trimmed_q30/1c_1_paired.fastq.gz \
#     ~/Downloads/macrogen_1c/trimmed_q30/1c_2_paired.fastq.gz \
#     ~/Downloads/macrogen_1c/filtering/contigs_final_clean.fasta

SAMPLE="$1"
R1="$2"
R2="$3"
ASSEMBLY="$4"

OUT_DIR="$HOME/Downloads/macrogen_1c/coverage/$SAMPLE"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

THREADS=$(sysctl -n hw.ncpu)
WINDOW_SIZE=1000   # bp per depth window -- adjust for finer/coarser resolution

echo "=== [$SAMPLE] Indexing assembly ==="
cp "$ASSEMBLY" ./assembly.fasta
bwa index assembly.fasta
samtools faidx assembly.fasta

echo ""
echo "=== [$SAMPLE] Aligning reads (bwa mem) and sorting ==="
bwa mem -t "$THREADS" assembly.fasta "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${SAMPLE}.sorted.bam" -
samtools index "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Alignment summary ==="
samtools flagstat "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Computing windowed depth (mosdepth, ${WINDOW_SIZE}bp windows) ==="
mosdepth --by "$WINDOW_SIZE" -t "$THREADS" "${SAMPLE}" "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Done ==="
echo "Depth windows: $OUT_DIR/${SAMPLE}.regions.bed.gz"
echo "Contig lengths: $OUT_DIR/assembly.fasta.fai"
