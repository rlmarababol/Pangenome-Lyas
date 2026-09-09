#!/bin/bash
set -euo pipefail

GENOME_DIR="$HOME/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/ani_genomes"
OUT_DIR="$HOME/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/paf_alignments"
mkdir -p "$OUT_DIR"

# Adjacent pairs in tree order: H, D, G, E, C, J, F, B, A, I
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
  echo "=== Aligning $g1 vs $g2 ==="
  minimap2 -x asm10 "$GENOME_DIR/${g1}.fasta" "$GENOME_DIR/${g2}.fasta" \
    > "$OUT_DIR/${g1}_vs_${g2}.paf" 2>/dev/null
  n_lines=$(wc -l < "$OUT_DIR/${g1}_vs_${g2}.paf")
  echo "  $n_lines alignment records"
done

echo ""
echo "=== All alignments complete ==="
ls -la "$OUT_DIR"
