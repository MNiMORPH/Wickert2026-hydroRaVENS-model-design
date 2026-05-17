#!/bin/bash
# Usage: bash archive_run.sh <run-name>
#
# Copies the configuration and outputs of the current Dakota run into
# runs/<run-name>/ for version-controlled storage.
#
# Files archived:
#   dakota.in, driver.py, params.yml, run_driver.sh  -- exact config
#   <config_template>  (resolved from params.yml)    -- hydroRaVENS config
#   evaluations.dat  (dakota.dat renamed)            -- all evaluations
#   dakota_log.txt   (dakota.out renamed)            -- Dakota log
#   best_fit.png     if present                      -- diagnostic plot

set -euo pipefail

NAME="${1:?Usage: bash archive_run.sh <run-name>}"
DEST="runs/${NAME}"

if [[ -d "$DEST" ]]; then
    echo "Error: $DEST already exists. Choose a different name." >&2
    exit 1
fi

mkdir -p "$DEST"

# Resolve config template from params.yml
CONFIG=$(python3 -c "
import yaml
with open('params.yml') as f:
    cfg = yaml.safe_load(f)
print(cfg['driver']['config_template'])
")

cp dakota.in              "$DEST/"
cp driver.py              "$DEST/"
cp params.yml             "$DEST/"
cp run_driver.sh          "$DEST/"
cp "$CONFIG"              "$DEST/"
cp dakota.dat             "$DEST/evaluations.dat"
cp dakota.out             "$DEST/dakota_log.txt"
[[ -f best_fit.png ]] && cp best_fit.png "$DEST/"

N=$(( $(wc -l < "$DEST/evaluations.dat") - 1 ))
echo "Archived to $DEST  ($N evaluations)"
