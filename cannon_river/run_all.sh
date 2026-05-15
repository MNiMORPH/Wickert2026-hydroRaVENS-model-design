#!/bin/bash
# Usage: bash run_all.sh <short-description>
# e.g.:  bash run_all.sh v1
#
# Runs model-selection experiments in series, each archived under
# runs/<timestamp>_<desc>/ inside its own experiment directory.
# Default: all M–S series. Override with --series M, N, O, P, Q, R, or S.

set -euo pipefail

DESC=""
SERIES="all"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --series) SERIES="$2"; shift 2 ;;
        *)        DESC="$1";   shift   ;;
    esac
done
DESC="${DESC:?Usage: bash run_all.sh <description> [--series M|N|O|P|Q|R|S|T|all]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

M_SERIES=( M00 M01 M02 M03 M04 M05 M06 )
N_SERIES=( N00 N01 N02 N03 N04 )
O_SERIES=( O00 O01 O02 O03 )
P_SERIES=( P00 P01 P02 P03 )
Q_SERIES=( Q00 Q01 Q02 Q03 )
R_SERIES=( R00 R01 R02 R03 R04 R05 )
S_SERIES=( S00 S01 )
T_SERIES=( T00 T01 )
U_SERIES=( U00 U01 )
case "$SERIES" in
    M)   EXPERIMENTS=( "${M_SERIES[@]}" ) ;;
    N)   EXPERIMENTS=( "${N_SERIES[@]}" ) ;;
    O)   EXPERIMENTS=( "${O_SERIES[@]}" ) ;;
    P)   EXPERIMENTS=( "${P_SERIES[@]}" ) ;;
    Q)   EXPERIMENTS=( "${Q_SERIES[@]}" ) ;;
    R)   EXPERIMENTS=( "${R_SERIES[@]}" ) ;;
    S)   EXPERIMENTS=( "${S_SERIES[@]}" ) ;;
    T)   EXPERIMENTS=( "${T_SERIES[@]}" ) ;;
    U)   EXPERIMENTS=( "${U_SERIES[@]}" ) ;;
    all) EXPERIMENTS=( "${M_SERIES[@]}" "${N_SERIES[@]}" "${O_SERIES[@]}" "${P_SERIES[@]}" "${Q_SERIES[@]}" "${R_SERIES[@]}" "${S_SERIES[@]}" "${T_SERIES[@]}" "${U_SERIES[@]}" ) ;;
    *)   echo "Unknown series '$SERIES'. Use M, N, O, P, Q, R, S, T, U, or all." >&2; exit 1 ;;
esac

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
echo "  collect AIC and goodness-of-fit across experiments."
echo "======================================================"
