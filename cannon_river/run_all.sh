#!/bin/bash
# Usage: bash run_all.sh <short-description> [--series LETTER] [--from EXP]
# e.g.:  bash run_all.sh v1
#        bash run_all.sh v1 --series V
#        bash run_all.sh v1 --series X --from X00F   # runs X00F, X01, X01F, ...
#
# Runs model-selection experiments in series, each archived under
# runs/<timestamp>_<desc>/ inside its own experiment directory.
# Default: all series. Override with --series M, N, O, P, Q, R, S, T, U, V, W, X, or Y.
# --from EXP skips experiments that sort before EXP (lexicographic; handles suffix
# letters such as X00F correctly).

set -euo pipefail

DESC=""
SERIES="all"
FROM=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --series) SERIES="$2"; shift 2 ;;
        --from)   FROM="$2";   shift 2 ;;
        *)        DESC="$1";   shift   ;;
    esac
done
DESC="${DESC:?Usage: bash run_all.sh <description> [--series M|N|O|P|Q|R|S|T|U|V|W|X|Y|all] [--from EXP]}"

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
V_SERIES=( V00 V01 V02 V03 V04 V05 V06 )
W_SERIES=( W00 W01 W02 )
X_SERIES=( X00 X00F X01 X01F X01Q X02 X02F X03 X03F )
Y_SERIES=( Y00 Y00F Y01 )
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
    V)   EXPERIMENTS=( "${V_SERIES[@]}" ) ;;
    W)   EXPERIMENTS=( "${W_SERIES[@]}" ) ;;
    X)   EXPERIMENTS=( "${X_SERIES[@]}" ) ;;
    Y)   EXPERIMENTS=( "${Y_SERIES[@]}" ) ;;
    all) EXPERIMENTS=( "${M_SERIES[@]}" "${N_SERIES[@]}" "${O_SERIES[@]}" "${P_SERIES[@]}" "${Q_SERIES[@]}" "${R_SERIES[@]}" "${S_SERIES[@]}" "${T_SERIES[@]}" "${U_SERIES[@]}" "${V_SERIES[@]}" "${W_SERIES[@]}" "${X_SERIES[@]}" "${Y_SERIES[@]}" ) ;;
    *)   echo "Unknown series '$SERIES'. Use M, N, O, P, Q, R, S, T, U, V, W, X, Y, or all." >&2; exit 1 ;;
esac

# Apply --from filter: include experiments at or after FROM (lexicographic order).
# Empty FROM means no filtering.
FILTERED=()
for EXP in "${EXPERIMENTS[@]}"; do
    [[ -z "$FROM" || "$EXP" > "$FROM" || "$EXP" == "$FROM" ]] && FILTERED+=("$EXP")
done
EXPERIMENTS=( "${FILTERED[@]}" )

if [[ ${#EXPERIMENTS[@]} -eq 0 ]]; then
    echo "No experiments match --series $SERIES --from $FROM." >&2
    exit 1
fi

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
