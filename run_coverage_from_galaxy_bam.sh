#!/bin/bash
set -euo pipefail

SAMPLE="LyasB"
GALAXY_BAM="/Users/ramjuliusmarababol/Downloads/Galaxy210-[Map with BWA on dataset 1, 3, and 162 (mapped reads in BAM format)].bam"

OUT_DIR="$HOME/Downloads/macrogen_1c/coverage/$SAMPLE"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

THREADS=$(sysctl -n hw.ncpu)
WINDOW_SIZE=1000

echo "=== [$SAMPLE] Checking BAM sort order ==="
SORT_ORDER=$(samtools view -H "$GALAXY_BAM" | grep "^@HD" | grep -o "SO:[a-zA-Z]*" || echo "SO:unknown")
echo "Detected: $SORT_ORDER"

if [[ "$SORT_ORDER" == "SO:coordinate" ]]; then
  echo "Already coordinate-sorted -- linking directly (no need to copy the large file)"
  ln -sf "$GALAXY_BAM" "${SAMPLE}.sorted.bam"
else
  echo "Not confirmed coordinate-sorted -- sorting now to be safe"
  samtools sort -@ "$THREADS" -o "${SAMPLE}.sorted.bam" "$GALAXY_BAM"
fi

echo ""
echo "=== [$SAMPLE] Indexing ==="
samtools index "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Alignment summary ==="
samtools flagstat "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Extracting contig lengths directly from BAM header ==="
samtools view -H "${SAMPLE}.sorted.bam" | grep "^@SQ" | \
  sed -E 's/.*SN:([^\t]+).*LN:([0-9]+).*/\1\t\2/' > assembly.fasta.fai
echo "Contigs found: $(wc -l < assembly.fasta.fai)"
head -5 assembly.fasta.fai

echo ""
echo "=== [$SAMPLE] Computing windowed depth (mosdepth, ${WINDOW_SIZE}bp windows) ==="
mosdepth --by "$WINDOW_SIZE" -t "$THREADS" "${SAMPLE}" "${SAMPLE}.sorted.bam"

echo ""
echo "=== [$SAMPLE] Done (using Galaxy-mapped BAM) ==="
echo "Depth windows: $OUT_DIR/${SAMPLE}.regions.bed.gz"
echo "Contig lengths: $OUT_DIR/assembly.fasta.fai"
