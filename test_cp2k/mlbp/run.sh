#!/usr/bin/env bash
#
# Run the CP2K pipeline on an ML-BP test case.
#
# Usage:
#   bash run.sh                # default: B3LYP
#   bash run.sh PBEh0125       # any of the seven functional tarballs
#
set -e

FUNC="${1:-B3LYP}"
TAG="$FUNC"

CIS_NTO_BIN="${CIS_NTO_BIN:-/home/hamid/cis_nto/build/bin}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

NGRID=60

# Resolve Molden / log paths by functional
case "$FUNC" in
    B3LYP)
        MOLDEN="b3lyp_opted-MOS-1_0.molden"
        LOG="out.log"
        ;;
    PBEh0125)
        MOLDEN="PBEh0125_50states/0125_PBE0.molden"
        LOG="PBEh0125_50states/out.log"
        ;;
    BHandHLYP)
        MOLDEN="BHandHLYP_50states/Bhand_opted-MOS-1_0.molden"
        LOG="BHandHLYP_50states/out.log"
        ;;
    PBE0)
        MOLDEN="PBE0_50_states/PBE0_OPTED-MOS-1_0.molden"
        LOG="PBE0_50_states/out.log"
        ;;
    CAM_B3LYP)
        MOLDEN="CAM_B3LYP_50_states/CAM_B3LYP_OPTED-MOS-1_0.molden"
        LOG="CAM_B3LYP_50_states/out.log"
        ;;
    wB97X)
        MOLDEN="wB97X_50_states/wB97x_opted-MOS-1_0.molden"
        LOG="wB97X_50_states/out.log"
        ;;
    wB97XD)
        MOLDEN="wb97XD_50_staes/w97x_d_opted-MOS-1_0.molden"
        LOG="wb97XD_50_staes/out.log"
        ;;
    *)
        # Generic fallback: find the first Molden under a matching directory
        dir=""
        for d in ${FUNC}*/ ; do
            if [ -d "$d" ]; then dir="$d"; break; fi
        done
        MOLDEN="$(find "$dir" -name "*_MOS-1_0.molden" -print -quit 2>/dev/null)"
        LOG="${dir}out.log"
        ;;
esac

[ -f "$MOLDEN" ] || { echo "ERROR: $MOLDEN missing. Run fetch.sh $FUNC first."; exit 1; }
[ -f "$LOG" ]    || { echo "ERROR: $LOG missing. Run fetch.sh $FUNC first.";    exit 1; }

mkdir -p cubes

echo "================================================================"
echo "  ML-BP $TAG test case"
echo "================================================================"
echo "  Molden: $MOLDEN"
echo "  Log:    $LOG"
echo

echo "  [1/5] convert.exe ..."
"$CIS_NTO_BIN/convert.exe" -in cp2k -out molden_cart \
    "${MOLDEN}:${LOG}" "std_${TAG}.molden" > /dev/null

echo "  [2/5] cis_nto.exe ..."
rm -f nto.molden
"$CIS_NTO_BIN/cis_nto.exe" -in cp2k "${MOLDEN}:${LOG}" \
    > "nto_${TAG}.log" 2>&1
[ -f nto.molden ] || { echo "ERROR: nto.molden not produced"; exit 1; }
mv nto.molden "nto_${TAG}.molden"

echo "  DONE: $TAG"
echo "  NTO log:    nto_${TAG}.log"
echo "  NTO Molden: nto_${TAG}.molden"

echo
echo "===== Reconstruction notice (from nto_${TAG}.log) ====="
grep -E "log references|molden contains|reconstruction method|dropping|FATAL" \
    "nto_${TAG}.log" || echo "(no reconstruction notice)"
echo
echo "Full reader output: nto_${TAG}.log"
echo "Structured record:  results.json (after 'python3 aggregate.py')"

echo
echo "===== Dominant NTO pair ====="
awk '
  /^ *Sym=/  { sym = $2 }
  /^ *Occup=/{ v = $2 + 0; if (v > 0.5) print sym, v }
' "nto_${TAG}.molden" | head -6
