#!/bin/bash
set -euo pipefail

GENOME_DIR="$HOME/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/ani_genomes"
OUT_DIR="$HOME/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/fastani_visual"
mkdir -p "$OUT_DIR"

declare -a PAIRS=(
  "LyasH LyasD"
  "LyasD LyasG"
  "LyasG LyasE"
  "LyasE LyasC"
  "LyasC LyasJ"
  "LyasJ LyasF"
  "LyasF LyasB"
  "LyasB LyasA"
  "LyasA LyasI"
)

for pair in "${PAIRS[@]}"; do
  read -r g1 g2 <<< "$pair"
  echo "=== $g1 (query) vs $g2 (reference) ==="
  fastANI -q "$GENOME_DIR/${g1}.fasta" -r "$GENOME_DIR/${g2}.fasta" \
    --visualize -o "$OUT_DIR/${g1}_vs_${g2}.out" 2>&1 | tail -3
  n_lines=$(wc -l < "$OUT_DIR/${g1}_vs_${g2}.out.visual")
  echo "  $n_lines mapping segments"
done

echo ""
echo "=== All fastANI visualize runs complete ==="
ls -la "$OUT_DIR"
