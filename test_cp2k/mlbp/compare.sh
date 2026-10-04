#!/usr/bin/env bash
#
# Verify ML-BP B3LYP output against reference values.
#
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

TAG="B3LYP"
EXPECTED_OCC="0.9940"
TOL="0.005"
MIN_NONZERO=20000

fail=0
echo "================================================================"
echo "  ML-BP B3LYP regression check"
echo "================================================================"
echo

echo "Checking dominant NTO occupation..."
nto="nto_${TAG}.molden"
[ -f "$nto" ] || { echo "  SKIP: $nto missing"; exit 1; }

occ=$(awk '/^ *Sym=/  { s = $2 }
           /^ *Occup=/{ v = $2 + 0; if (s ~ /^2/ && v > best) best = v }
           END{ printf "%.4f", best }' "$nto")
if awk -v a="$occ" -v b="$EXPECTED_OCC" -v t="$TOL" \
    'BEGIN{ d = a - b; exit !(d < t && d > -t) }'; then
    echo "  PASS  occ=$occ   (expected $EXPECTED_OCC ± $TOL)"
else
    echo "  FAIL  occ=$occ   (expected $EXPECTED_OCC ± $TOL)"
    fail=$((fail + 1))
fi

echo
echo "Checking cube integrals and nonzero counts..."
for orb in s1_occ s1_virt homo lumo; do
    cube="cubes/${orb}_${TAG}.cube"
    [ -f "$cube" ] || { echo "  SKIP: $cube"; continue; }

    read -r integral nonzero < <(python3 - "$cube" <<'PYEOF'
import sys, numpy as np
lines = open(sys.argv[1]).readlines()
natoms = int(lines[2].split()[0])
dx = float(lines[3].split()[1]); dy = float(lines[4].split()[2]); dz = float(lines[5].split()[3])
vals = []
for line in lines[6+natoms:]:
    for tok in line.split():
        try: vals.append(float(tok))
        except ValueError: pass
arr = np.array(vals)
print(f"{float(np.sum(arr*arr))*dx*dy*dz:.6f} {int(np.sum(np.abs(arr) > 1e-6))}")
PYEOF
)

    ok=1
    awk -v a="$integral" 'BEGIN{ d = a - 1.0; exit !(d < 1e-3 && d > -1e-3) }' || ok=0
    [ "$nonzero" -ge "$MIN_NONZERO" ] || ok=0

    if [ "$ok" -eq 1 ]; then
        echo "  PASS  $orb   ∫|ψ|²=$integral   nonzero=$nonzero"
    else
        echo "  FAIL  $orb   ∫|ψ|²=$integral   nonzero=$nonzero"
        fail=$((fail + 1))
    fi
done

echo
if [ "$fail" -eq 0 ]; then echo "  ALL CHECKS PASSED"; exit 0
else echo "  $fail CHECK(S) FAILED"; exit 1; fi
