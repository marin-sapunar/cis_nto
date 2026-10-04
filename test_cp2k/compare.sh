#!/usr/bin/env bash
#
# Regression check for the CP2K adamantane-family test case.
#
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

declare -A EXPECTED_OCC=(
    ["adamantane"]="0.9976"
    ["1methyladamantane"]="0.9980"
    ["13dimethyladamantane"]="0.9961"
    ["135trimethyladamantane"]="0.9976"
)
TAGS=(adamantane 1methyladamantane 13dimethyladamantane 135trimethyladamantane)
MIN_NONZERO=2000

fail=0
echo "================================================================"
echo "  CP2K adamantane-family regression check"
echo "================================================================"
echo

echo "Checking dominant S1 NTO occupations..."
for tag in "${TAGS[@]}"; do
    nto="nto_${tag}.molden"
    [ -f "$nto" ] || { echo "  SKIP: $nto not found"; continue; }
    occ=$(awk '/^ *Sym=/  { s = $2 }
               /^ *Occup=/{ v = $2 + 0; if (s ~ /^2/ && v > best) best = v }
               END{ printf "%.4f", best }' "$nto")
    expected="${EXPECTED_OCC[$tag]}"
    if awk -v a="$occ" -v b="$expected" \
        'BEGIN{ d = a - b; exit !(d < 1e-3 && d > -1e-3) }'; then
        echo "  PASS  $tag   occ=$occ   (expected $expected)"
    else
        echo "  FAIL  $tag   occ=$occ   (expected $expected)"
        fail=$((fail + 1))
    fi
done

echo
echo "Checking cube integrals and nonzero counts..."
for tag in "${TAGS[@]}"; do
    for orb in s1_occ s1_virt homo lumo; do
        cube="cubes/${orb}_${tag}.cube"
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
            echo "  PASS  $orb $tag   ∫|ψ|²=$integral   nonzero=$nonzero"
        else
            echo "  FAIL  $orb $tag   ∫|ψ|²=$integral   nonzero=$nonzero"
            fail=$((fail + 1))
        fi
    done
done

echo
if [ "$fail" -eq 0 ]; then echo "  ALL CHECKS PASSED"; exit 0
else echo "  $fail CHECK(S) FAILED"; exit 1; fi
