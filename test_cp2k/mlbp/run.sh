#!/usr/bin/env bash
#
# Run the CP2K pipeline on the ML-BP B3LYP test case.
#
# Requires: fetch.sh has been run (Molden + log present).
#
set -e

CIS_NTO_BIN="${CIS_NTO_BIN:-/home/hamid/cis_nto/build/bin}"
HERE="$(cd "$(dirname "$0")" && pwd)"
TEST_ROOT="$(cd "$HERE/.." && pwd)"

cd "$HERE"

MOLDEN="b3lyp_opted-MOS-1_0.molden"
LOG="out.log"
TAG="B3LYP"

if [ ! -f "$MOLDEN" ]; then
    echo "ERROR: $MOLDEN not found. Run fetch.sh first."
    exit 1
fi
if [ ! -f "$LOG" ]; then
    echo "ERROR: $LOG not found. Run fetch.sh first."
    exit 1
fi

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

if [ ! -f nto.molden ]; then
    echo "  ERROR: nto.molden not produced"
    exit 1
fi
mv nto.molden "nto_${TAG}.molden"

echo "  [3/5] extracting S1 NTO pair ..."
dom_sym=$(awk '
    /Sym=/  { sym = $2 }
    /Occup/ { v = $2 + 0; if (sym ~ /^2/ && v > 0.5) { print sym; exit } }
' "nto_${TAG}.molden")
if [ -z "$dom_sym" ]; then
    echo "  ERROR: no dominant occupied NTO found"
    exit 1
fi

virt_sym=$(awk '
    /Sym=/  { sym = $2 }
    /Occup/ { v = $2 + 0; if (sym ~ /^1/ && v > 0.5) { print sym; exit } }
' "nto_${TAG}.molden")
if [ -z "$virt_sym" ]; then
    echo "  ERROR: no dominant virtual NTO found"
    exit 1
fi

echo "    dominant occupied NTO: $dom_sym"
echo "    dominant virtual  NTO: $virt_sym"

python3 "$TEST_ROOT/extract_nto.py" "nto_${TAG}.molden" "s1_occ_${TAG}.molden"  "$dom_sym"
python3 "$TEST_ROOT/extract_nto.py" "nto_${TAG}.molden" "s1_virt_${TAG}.molden" "$virt_sym"

echo "  [4/5] extracting HOMO/LUMO ..."
python3 - "$TAG" <<'PYEOF'
import sys
tag = sys.argv[1]
lines = open(f'std_{tag}.molden').readlines()
syms = [l.split()[1] for l in lines if l.strip().startswith('Sym=')]
open(f'homo_sym_{tag}.txt','w').write(syms[159] if len(syms) > 159 else '')
open(f'lumo_sym_{tag}.txt','w').write(syms[160] if len(syms) > 160 else '')
PYEOF

homo_sym=$(cat "homo_sym_${TAG}.txt")
lumo_sym=$(cat "lumo_sym_${TAG}.txt")
python3 "$TEST_ROOT/extract_nto.py" "std_${TAG}.molden" "homo_${TAG}.molden" "$homo_sym"
python3 "$TEST_ROOT/extract_nto.py" "std_${TAG}.molden" "lumo_${TAG}.molden" "$lumo_sym"

echo "  [5/5] generating cube files ..."
for orb in homo lumo s1_occ s1_virt; do
    python3 "$TEST_ROOT/molden_to_cube.py" \
        "${orb}_${TAG}.molden" "cubes/${orb}_${TAG}.cube" 60
done

echo "  DONE: $TAG"
echo
echo "Cube files are in cubes/"
ls -la cubes/