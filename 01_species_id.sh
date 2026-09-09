#!/bin/bash
# 01_species_id.sh
# Species-level identity confirmation via average nucleotide identity
# against the L. yasudae type strain F1(T) -- Supplementary Methods S8.
# Usage: ./01_species_id.sh <query_assembly.fasta> <F1T_reference.fasta> <output.ani>

set -euo pipefail

QUERY="$1"
REFERENCE="$2"
OUTPUT="$3"

echo "=== Species identification (FastANI vs. L. yasudae F1T) ==="
fastANI -q "$QUERY" -r "$REFERENCE" -o "$OUTPUT"

echo "Done. Result: $OUTPUT"
cat "$OUTPUT"
