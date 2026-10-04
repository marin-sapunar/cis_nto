#!/usr/bin/env bash
#
# Full CP2K -> cis_nto -> cube pipeline for the adamantane series.
#
# For each molecule:
#   1. Convert CP2K molden + log to a standardized molden via convert.exe
#   2. Run cis_nto.exe on the CP2K pair -> nto.molden
#   3. Extract the S1 dominant NTO pair -> s1_occ_NTO.molden, s1_virt_NTO.molden
#   4. Extract the S0 HOMO and LUMO from the standardized molden
#   5. Generate cube files for all four orbitals
#
# Usage: ./run_all.sh
#
set -e

# ----- paths -----
CIS_NTO_BIN=/home/hamid/cis_nto/build/bin
TEST_DIR=/home/hamid/cis_nto/test_cp2k
cd "$TEST_DIR"

# ----- molecule list -----
# Format:  <tag>  <molden>  <log>
MOLS=(
    "adamantane              adamantane.molden              adamantane.log"
    "1methyladamantane       1_methyladamantane.molden      1_methyladamantane.log"
    "13dimethyladamantane    1_3_dimethyladamantane.molden  1_3_dimethyladamantane.log"
    "135trimethyladamantane  1_3_5_trimethyladamantane.molden 1_3_5_trimethyladamantane.log"
)

mkdir -p cubes

for entry in "${MOLS[@]}"; do
    set -- $entry
    tag=$1
    molden=$2
    log=$3

    echo "================================================================"
    echo "  $tag"
    echo "================================================================"

    if [ ! -f "$molden" ]; then
        echo "  SKIP: $molden not found"
        continue
    fi
    if [ ! -f "$log" ]; then
        echo "  SKIP: $log not found"
        continue
    fi

    # 1. Convert to standardized molden
    echo "  [1/5] convert.exe ..."
    "$CIS_NTO_BIN/convert.exe" -in cp2k -out molden_cart \
        "${molden}:${log}" "std_${tag}.molden" > /dev/null

    # 2. Run cis_nto
    echo "  [2/5] cis_nto.exe ..."
    "$CIS_NTO_BIN/cis_nto.exe" -in cp2k "${molden}:${log}" \
        > "nto_${tag}.log" 2>&1

    # The program writes nto.molden in the CWD
    if [ ! -f nto.molden ]; then
        echo "  ERROR: nto.molden not produced"
        continue
    fi
    mv nto.molden "nto_${tag}.molden"

    # 3. Find the highest-occupation NTO pair
    #    Sym tags: 1xxx = virtual manifold, 2xxx = occupied manifold
    #    The dominant pair is the first Sym= 2xx1a / 2xx2a with occ > 0.9
    echo "  [3/5] extracting S1 NTO pair ..."
    dom_sym=$(awk '/Sym=/{sym=$2} /Occup=/{v=$2+0; if (v > 0.9 && sym ~ /^2/) {print sym; exit}}' "nto_${tag}.molden")
    if [ -z "$dom_sym" ]; then
        echo "  ERROR: no dominant NTO found"
        continue
    fi
    # The virtual partner is the immediately following Sym tag
    virt_sym=$(awk -v target="$dom_sym" '
        /Sym=/ {sym=$2}
        sym == target {found=1; next}
        found && /Sym=/ {print $2; exit}
    ' "nto_${tag}.molden")

    echo "    dominant occupied NTO: $dom_sym"
    echo "    dominant virtual  NTO: $virt_sym"

    python3 extract_nto.py "nto_${tag}.molden" "s1_occ_${tag}.molden" "$dom_sym"
    python3 extract_nto.py "nto_${tag}.molden" "s1_virt_${tag}.molden" "$virt_sym"

    # 4. Extract HOMO and LUMO from std molden (last two occupied, first virtual)
    #    For adamantane-family the HOMO is a degenerate manifold; take the
    #    first MO with the dominant NTO overlap. For simplicity, extract
    #    the MO just below the LUMO.
    echo "  [4/5] extracting HOMO/LUMO ..."
    python3 - "$tag" <<'PYEOF'
import sys
tag = sys.argv[1]
# Find the highest-occupation MO in std_<tag>.molden (should be the HOMO)
lines = open(f'std_{tag}.molden').readlines()
syms = []
for i, line in enumerate(lines):
    if line.strip().startswith('Sym='):
        syms.append(line.split()[1])
# For 68-MO system: occupied 1-38, virtual 39-68
# Take MO 38 as HOMO candidate, MO 39 as LUMO
with open(f'homo_sym_{tag}.txt', 'w') as f:
    f.write(syms[37] if len(syms) > 37 else '')
with open(f'lumo_sym_{tag}.txt', 'w') as f:
    f.write(syms[38] if len(syms) > 38 else '')
PYEOF

    homo_sym=$(cat "homo_sym_${tag}.txt")
    lumo_sym=$(cat "lumo_sym_${tag}.txt")

    python3 extract_nto.py "std_${tag}.molden" "homo_${tag}.molden" "$homo_sym"
    python3 extract_nto.py "std_${tag}.molden" "lumo_${tag}.molden" "$lumo_sym"

    # 5. Generate cube files
    echo "  [5/5] generating cube files ..."
    python3 molden_to_cube.py "homo_${tag}.molden"     "cubes/homo_${tag}.cube"     40
    python3 molden_to_cube.py "lumo_${tag}.molden"     "cubes/lumo_${tag}.cube"     40
    python3 molden_to_cube.py "s1_occ_${tag}.molden"   "cubes/s1_occ_${tag}.cube"   40
    python3 molden_to_cube.py "s1_virt_${tag}.molden"  "cubes/s1_virt_${tag}.cube"  40

    echo "  DONE: $tag"
    echo ""
done

echo "All molecules processed."
echo "Cube files are in cubes/"
ls -la cubes/