#!/usr/bin/env bash
#
# Regression check for the CP2K test case.
#
# Verifies that:
#   1. Each molecule's nto_<tag>.molden has the expected dominant NTO pair.
#   2. Each generated cube satisfies ∫|ψ|² dV = 1.
#
# Run from this directory after ./run_all.sh has completed.
#
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

# Expected dominant NTO occupations for S1 of each molecule.
# These were measured from the reference run; they should be reproduced
# to within 1e-3 on any correct execution of the pipeline.
declare -A EXPECTED=(
    ["adamantane"]="0.9976"
    ["1methyladamantane"]="0.9980"
    ["13dimethyladamantane"]="0.9961"
    ["135trimethyladamantane"]="0.9976"
)

TAGS=(adamantane 1methyladamantane 13dimethyladamantane 135trimethyladamantane)

echo "================================================================"
echo "  CP2K test regression check"
echo "================================================================"
echo

fail=0

# ---- Check 1: NTO occupations ----------------------------------------------
echo "Checking dominant NTO occupations..."
for tag in "${TAGS[@]}"; do
    nto="nto_${tag}.molden"
    if [ ! -f "$nto" ]; then
        echo "  SKIP: $nto not found (run ./run_all.sh first)"
        continue
    fi

    # Highest occupation in the 2xxx manifold (occupied NTO)
    occ=$(awk '/Sym=/{s=$2} /Occup=/{v=$2+0; if (s ~ /^2/ && v > best) { best = v } } END{printf "%.4f", best}' "$nto")

    expected="${EXPECTED[$tag]}"
    delta=$(awk -v a="$occ" -v b="$expected" 'BEGIN{printf "%.4f", a-b}')

    if [ "$(awk -v d="$delta" 'BEGIN{exit !(d<0.001 && d>-0.001)}' ; echo $?)" = "0" ]; then
        echo "  PASS  $tag   occ=$occ   (expected $expected)"
    else
        echo "  FAIL  $tag   occ=$occ   (expected $expected, delta $delta)"
        fail=$((fail + 1))
    fi
done

echo

# ---- Check 2: Cube integrals ----------------------------------------------
echo "Checking cube L2 normalization..."
for tag in "${TAGS[@]}"; do
    for orb in s1_occ s1_virt homo lumo; do
        cube="cubes/${orb}_${tag}.cube"
        if [ ! -f "$cube" ]; then
            echo "  SKIP: $cube not found"
            continue
        fi

        integral=$(python3 - "$cube" <<'PYEOF'
import sys
import numpy as np

fn = sys.argv[1]
lines = open(fn).readlines()
natoms = int(lines[2].split()[0])
nx = int(lines[3].split()[0]); dx = float(lines[3].split()[1])
ny = int(lines[4].split()[0]); dy = float(lines[4].split()[2])
nz = int(lines[5].split()[0]); dz = float(lines[5].split()[3])
vals = []
for line in lines[6+natoms:]:
    for tok in line.split():
        try: vals.append(float(tok))
        except ValueError: pass
dV = dx * dy * dz
integral = float(np.sum(np.array(vals)**2)) * dV
print(f"{integral:.6f}")
PYEOF
)

        delta=$(awk -v a="$integral" 'BEGIN{printf "%.6f", a-1.0}')
        if [ "$(awk -v d="$delta" 'BEGIN{exit !(d<0.001 && d>-0.001)}' ; echo $?)" = "0" ]; then
            echo "  PASS  $orb $tag   ∫|ψ|² dV = $integral"
        else
            echo "  FAIL  $orb $tag   ∫|ψ|² dV = $integral  (expected 1.0)"
            fail=$((fail + 1))
        fi
    done
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