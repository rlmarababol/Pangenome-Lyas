#!/bin/bash
# 05_scaffold.sh
# Reference-guided scaffolding (RagTag) -- Supplementary Methods S8.
# Usage: ./05_scaffold.sh <reference.fasta> <contigs.fasta> <output_dir>

set -euo pipefail

REFERENCE="$1"
CONTIGS="$2"
OUTDIR="$3"

echo "=== Reference-guided scaffolding (RagTag) ==="
ragtag.py scaffold "$REFERENCE" "$CONTIGS" -o "$OUTDIR"

echo "Done. Scaffolded assembly: $OUTDIR/ragtag.scaffold.fasta"
