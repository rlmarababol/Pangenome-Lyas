#!/bin/bash
# 03_mlst.sh
# Multilocus sequence typing against the Leptospira PubMLST scheme --
# Supplementary Methods S8.
# Usage: ./03_mlst.sh <assembly.fasta> <output.tsv>

set -euo pipefail

ASSEMBLY="$1"
OUTPUT="$2"

echo "=== MLST (Leptospira scheme) ==="
mlst --scheme leptospira "$ASSEMBLY" > "$OUTPUT"

echo "Done. Result: $OUTPUT"
cat "$OUTPUT"
echo ""
echo "NOTE: allele calls prefixed with '~' indicate approximate (inexact)"
echo "matches to the nearest catalogued allele, not exact matches."
