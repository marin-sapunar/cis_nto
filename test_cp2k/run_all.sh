#!/usr/bin/env bash
#
# Full CP2K -> cis_nto -> cube pipeline for the adamantane series.
#
set -e

CIS_NTO_BIN="${CIS_NTO_BIN:-/home/hamid/cis_nto/build/bin}"
TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$TEST_DIR"

NGRID=40

MOLS=(
    "adamantane              adamantane.molden              adamantane.log"
    "1methyladamantane       1_methyladamantane.molden      1_methyladamantane.log"
    "13dimethyladamantane    1_3_dimethyladamantane.molden  1_3_dimethyladamantane.log"
    "135trimethyladamantane  1_3_5_trimethyladamantane.molden 1_3_5_trimethyladamantane.log"
)

mkdir -p cubes

for entry in "${MOLS[@]}"; do
    set -- $entry
    tag=$1; molden=$2; log=$3

    echo "================================================================"
    echo "  $tag"
    echo "================================================================"

    [ -f "$molden" ] || { echo "  SKIP: $molden not found"; continue; }
    [ -f "$log" ]    || { echo "  SKIP: $log not found";    continue; }

    echo "  [1/5] convert.exe ..."
    "$CIS_NTO_BIN/convert.exe" -in cp2k -out molden_cart \
        "${molden}:${log}" "std_${tag}.molden" > /dev/null

    echo "  [2/5] cis_nto.exe ..."
    "$CIS_NTO_BIN/cis_nto.exe" -in cp2k "${molden}:${log}" \
        > "nto_${tag}.log" 2>&1
    [ -f nto.molden ] || { echo "  ERROR: nto.molden not produced"; continue; }
    mv nto.molden "nto_${tag}.molden"

    echo "  [3/5] extracting S1 NTO pair ..."
    dom_sym=$(awk '
        /^ *Sym=/  { sym = $2 }
        /^ *Occup=/{ v = $2 + 0; if (sym ~ /^2/ && v > 0.5) { print sym; exit } }
    ' "nto_${tag}.molden")
    virt_sym=$(awk '
        /^ *Sym=/  { sym = $2 }
        /^ *Occup=/{ v = $2 + 0; if (sym ~ /^1/ && v > 0.5) { print sym; exit } }
    ' "nto_${tag}.molden")

    [ -n "$dom_sym" ]  || { echo "  ERROR: no dominant occupied NTO"; continue; }
    [ -n "$virt_sym" ] || { echo "  ERROR: no dominant virtual NTO";  continue; }

    echo "    occupied: $dom_sym   virtual: $virt_sym"

    python3 extract_nto.py "nto_${tag}.molden" "s1_occ_${tag}.molden"  "$dom_sym"
    python3 extract_nto.py "nto_${tag}.molden" "s1_virt_${tag}.molden" "$virt_sym"

    echo "  [4/5] extracting HOMO/LUMO by index ..."
    python3 extract_nto.py "std_${tag}.molden" "homo_${tag}.molden" 38
    python3 extract_nto.py "std_${tag}.molden" "lumo_${tag}.molden" 39

    echo "  [5/5] generating cube files ..."
    for orb in homo lumo s1_occ s1_virt; do
        python3 molden_to_cube.py "${orb}_${tag}.molden" \
            "cubes/${orb}_${tag}.cube" "$NGRID"
    done

    echo "  DONE: $tag"
    echo ""
done

echo "All molecules processed."
ls -la cubes/
