#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${1:-$ROOT_DIR/tools/examples/universe_state.json}"
OUTPUT="${2:-$ROOT_DIR/out/content}"
SEED="${3:-123}"

python -m tools.stellar_content_pipeline --input "$INPUT" --output "$OUTPUT" --seed "$SEED"

echo "Done. See: $OUTPUT"
