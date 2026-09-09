#!/bin/bash
set -euo pipefail

GFF3_DIR="/Users/ramjuliusmarababol/Downloads/GFF3"
OUT_DIR="$HOME/Downloads/macrogen_1c/pangenome/pangenome_results/synteny"
GENOME_DIR="$OUT_DIR/ani_genomes"

mkdir -p "$GENOME_DIR"
cd "$GFF3_DIR"

echo "=== Extracting genome FASTAs from GFF3 (fresh, to be safe) ==="
for code in LyasA LyasB LyasC LyasD LyasE LyasF LyasG LyasH LyasI LyasJ; do
  awk '/^##FASTA/{found=1; next} found' "${code}.gff3" > "$GENOME_DIR/${code}.fasta"
done

ls -la "$GENOME_DIR"

cd "$OUT_DIR"
ls "$GENOME_DIR"/*.fasta > genome_list.txt

echo ""
echo "=== Running fastANI (all-vs-all) ==="
fastANI --ql genome_list.txt --rl genome_list.txt \
  -o ani_all_vs_all.tsv -t "$(sysctl -n hw.ncpu)"

echo ""
echo "=== Raw fastANI output (head) ==="
head -10 ani_all_vs_all.tsv
wc -l ani_all_vs_all.tsv

echo ""
echo "Saved: $OUT_DIR/ani_all_vs_all.tsv"
