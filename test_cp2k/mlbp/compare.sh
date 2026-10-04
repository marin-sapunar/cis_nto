#!/usr/bin/env bash
#
# Verify ML-BP B3LYP output against reference values.
#
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

TAG="B3LYP"
EXPECTED_OCC="0.9940"   # dominant S1 NTO occupation for B3LYP
TOL="0.005"

fail=0

echo "================================================================"
echo "  ML-BP B3LYP regression check"
echo "================================================================"
echo

# ---- Check 1: dominant NTO occupation --------------------------------------
echo "Checking dominant NTO occupation..."
nto="nto_${TAG}.molden"
if [ ! -f "$nto" ]; then
    echo "  SKIP: $nto not found (run ./run.sh first)"
    exit 1
fi

occ=$(awk '/Sym=/{s=$2} /Occup=/{v=$2+0; if (s ~ /^2/ && v > best) { best = v } } END{printf "%.4f", best}' "$nto")
delta=$(awk -v a="$occ" -v b="$EXPECTED_OCC" 'BEGIN{printf "%.4f", a-b}')
if [ "$(awk -v d="$delta" -v t="$TOL" 'BEGIN{exit !(d<t && d>-t)}' ; echo $?)" = "0" ]; then
    echo "  PASS  occ=$occ   (expected $EXPECTED_OCC ± $TOL)"
else
    echo "  FAIL  occ=$occ   (expected $EXPECTED_OCC ± $TOL, delta $delta)"
    fail=$((fail + 1))
fi

echo

# ---- Check 2: cube L2 normalization ----------------------------------------
echo "Checking cube L2 normalization..."
for orb in s1_occ s1_virt homo lumo; do
    cube="cubes/${orb}_${TAG}.cube"
    if [ ! -f "$cube" ]; then
        echo "  SKIP: $cube not found"
        continue
    fi

    integral=$(python3 - "$cube" <<'PYEOF'
import sys, numpy as np
lines = open(sys.argv[1]).readlines()
natoms = int(lines[2].split()[0])
dx = float(lines[3].split()[1]); dy = float(lines[4].split()[2]); dz = float(lines[5].split()[3])
vals = []
for line in lines[6+natoms:]:
    for tok in line.split():
        try: vals.append(float(tok))
        except ValueError: pass
print(f"{float(np.sum(np.array(vals)**2)) * dx*dy*dz:.6f}")
PYEOF
)

    delta=$(awk -v a="$integral" 'BEGIN{printf "%.6f", a-1.0}')
    if [ "$(awk -v d="$delta" 'BEGIN{exit !(d<0.001 && d>-0.001)}' ; echo $?)" = "0" ]; then
        echo "  PASS  $orb   ∫|ψ|² dV = $integral"
    else
        echo "  FAIL  $orb   ∫|ψ|² dV = $integral  (expected 1.0)"
        fail=$((fail + 1))
    fi
done

echo
if [ "$fail" -eq 0 ]; then
    echo "================================================================"
    echo "  ALL CHECKS PASSED"
    echo "================================================================"
    exit 0
else
    echo "================================================================"
    echo "  $fail CHECK(S) FAILED"
    echo "================================================================"
    exit 1
fi