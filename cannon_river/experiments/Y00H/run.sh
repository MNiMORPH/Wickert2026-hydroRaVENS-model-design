#!/bin/bash
# Usage: bash run.sh <short-description>
# e.g.:  bash run.sh kge_2res_hmax
#
# Cleans previous ephemeral outputs, runs Dakota, generates the best-fit
# diagnostic plot, and archives all results to runs/<timestamp>_<desc>/.
# The timestamp prefix guarantees no run overwrites another.

set -euo pipefail

DESC="${1:?Usage: bash run.sh <short-description>  e.g. kge_2res_hmax}"
TIMESTAMP=$(date +%Y-%m-%d_%H%M%S)
RUN_NAME="${TIMESTAMP}_${DESC}"

DAKOTA=/home/awickert/anaconda3/envs/dakota-env/bin/dakota
PYTHON=/home/awickert/anaconda3/envs/dakota-env/bin/python

echo "=== Run: $RUN_NAME ==="

# Regenerate dakota.in from params.yml to keep them in sync
$PYTHON generate_dakota_in.py

# Clean previous ephemeral outputs
rm -rf out dakota.dat dakota.out dakota.rst fort.13 LHS_*.out

# Optimise
$DAKOTA -i dakota.in -o dakota.out

# Save figure without showing it so we can archive before blocking on display.
if $PYTHON plot_best.py --save best_fit.png --no-show; then
    echo "Best-fit plot saved."
else
    echo "Warning: plot_best.py failed; archiving without plot." >&2
fi

# Archive while dakota.dat / best_fit.png still belong to this run
bash archive_run.sh "$RUN_NAME"

echo "=== Archived to runs/$RUN_NAME ==="

# Open the archived figure non-blocking so new runs are unaffected
[[ -f "runs/$RUN_NAME/best_fit.png" ]] && xdg-open "runs/$RUN_NAME/best_fit.png" &
