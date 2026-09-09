#!/bin/bash
# 02_strain_ani.sh
# Strain-level relatedness between two isolates via pairwise ANI --
# Supplementary Methods S8. Used to compare LyasA vs. LyasB directly.
# Usage: ./02_strain_ani.sh <isolate1.fasta> <isolate2.fasta> <output.ani>

set -euo pipefail

ISOLATE1="$1"
ISOLATE2="$2"
OUTPUT="$3"

echo "=== Strain-level ANI comparison ==="
fastANI -q "$ISOLATE1" -r "$ISOLATE2" -o "$OUTPUT"

echo "Done. Result: $OUTPUT"
cat "$OUTPUT"
