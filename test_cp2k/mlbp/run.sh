#!/usr/bin/env bash
#
# Run the CP2K pipeline on the ML-BP B3LYP test case.
#
set -e

CIS_NTO_BIN="${CIS_NTO_BIN:-/home/hamid/cis_nto/build/bin}"
HERE="$(cd "$(dirname "$0")" && pwd)"
TEST_ROOT="$(cd "$HERE/.." && pwd)"
cd "$HERE"

NGRID=60
MOLDEN="b3lyp_opted-MOS-1_0.molden"
LOG="out.log"
TAG="B3LYP"

[ -f "$MOLDEN" ] || { echo "ERROR: $MOLDEN missing. Run fetch.sh."; exit 1; }
[ -f "$LOG" ]    || { echo "ERROR: $LOG missing. Run fetch.sh.";    exit 1; }

mkdir -p cubes

echo "================================================================"
echo "  ML-BP B3LYP test case"
echo "================================================================"

echo "  [1/5] convert.exe ..."
"$CIS_NTO_BIN/convert.exe" -in cp2k -out molden_cart \
    "${MOLDEN}:${LOG}" "std_${TAG}.molden" > /dev/null

echo "  [2/5] cis_nto.exe ..."
"$CIS_NTO_BIN/cis_nto.exe" -in cp2k "${MOLDEN}:${LOG}" \
    > "nto_${TAG}.log" 2>&1
[ -f nto.molden ] || { echo "ERROR: nto.molden not produced"; exit 1; }
mv nto.molden "nto_${TAG}.molden"

echo "  [3/5] extracting S1 NTO pair ..."
dom_sym=$(awk '
    /^ *Sym=/  { sym = $2 }
    /^ *Occup=/{ v = $2 + 0; if (sym ~ /^2/ && v > 0.5) { print sym; exit } }
' "nto_${TAG}.molden")
virt_sym=$(awk '
    /^ *Sym=/  { sym = $2 }
    /^ *Occup=/{ v = $2 + 0; if (sym ~ /^1/ && v > 0.5) { print sym; exit } }
' "nto_${TAG}.molden")

[ -n "$dom_sym" ]  || { echo "ERROR: no dominant occupied NTO"; exit 1; }
[ -n "$virt_sym" ] || { echo "ERROR: no dominant virtual NTO";  exit 1; }

echo "    occupied: $dom_sym   virtual: $virt_sym"

python3 "$TEST_ROOT/extract_nto.py" "nto_${TAG}.molden" "s1_occ_${TAG}.molden"  "$dom_sym"
python3 "$TEST_ROOT/extract_nto.py" "nto_${TAG}.molden" "s1_virt_${TAG}.molden" "$virt_sym"

echo "  [4/5] extracting HOMO/LUMO by index (160, 161) ..."
python3 "$TEST_ROOT/extract_nto.py" "std_${TAG}.molden" "homo_${TAG}.molden" 160
python3 "$TEST_ROOT/extract_nto.py" "std_${TAG}.molden" "lumo_${TAG}.molden" 161

echo "  [5/5] generating cube files ..."
for orb in homo lumo s1_occ s1_virt; do
    python3 "$TEST_ROOT/molden_to_cube.py" \
        "${orb}_${TAG}.molden" "cubes/${orb}_${TAG}.cube" "$NGRID"
done

echo "  DONE: $TAG"
ls -la cubes/
