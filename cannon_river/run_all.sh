#!/bin/bash
# Usage: bash run_all.sh <short-description>
# e.g.:  bash run_all.sh v1
#
# Runs all model-selection experiments (M00–M06) in series, each archived
# under runs/<timestamp>_<desc>/ inside its own experiment directory.

set -euo pipefail

DESC="${1:?Usage: bash run_all.sh <short-description>  e.g. v1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENTS=( M00 M01 M02 M03 M04 M05 M06 )

for EXP in "${EXPERIMENTS[@]}"; do
    echo ""
    echo "======================================================"
    echo "  Starting $EXP  ($(date '+%Y-%m-%d %H:%M:%S'))"
    echo "======================================================"
    cd "$SCRIPT_DIR/experiments/$EXP"
    bash run.sh "$DESC"
done

echo ""
echo "======================================================"
echo "  All experiments complete  ($(date '+%Y-%m-%d %H:%M:%S'))"
echo "======================================================"
echo "  Run 'python compare_aic.py' from cannon_river/ to"
echo "  collect AIC and goodness-of-fit across M00–M06."
echo "======================================================"
