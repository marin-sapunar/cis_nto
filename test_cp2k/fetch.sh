#!/usr/bin/env bash
#
# Fetch the adamantane-family test inputs from the AkimovLab repository.
#
# Source: https://github.com/AkimovLab/Project_adamantanes
# Path:   Static_Optimization_Structure/
#
set -eu

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

BASE="https://raw.githubusercontent.com/AkimovLab/Project_adamantanes/main/Static_Optimization_Structure"

# "<upstream_path> <local_name>" pairs
PAIRS="
adamantane/0_molden.molden adamantane.molden
adamantane/0_output.log adamantane.log
1_methyladamantane/1_molden.molden 1_methyladamantane.molden
1_methyladamantane/1_output.log 1_methyladamantane.log
1_3_dimethyladamantane/2_molden.molden 1_3_dimethyladamantane.molden
1_3_dimethyladamantane/2_output.log 1_3_dimethyladamantane.log
1_3_5_trimethyladamantane/3_molden.molden 1_3_5_trimethyladamantane.molden
1_3_5_trimethyladamantane/3_output.log 1_3_5_trimethyladamantane.log
"

echo "Fetching adamantane-family test inputs..."
echo

echo "$PAIRS" | while read -r upstream local; do
    [ -z "$upstream" ] && continue

    if [ -f "$local" ]; then
        echo "  present: $local"
        continue
    fi

    url="$BASE/$upstream"
    if curl -fsSL -o "$local" "$url"; then
        echo "  fetched: $local"
    else
        rm -f "$local"
        echo "  MISSING: $local"
    fi
done

# Verify all 8 files exist
missing=0
for f in adamantane.molden adamantane.log \
         1_methyladamantane.molden 1_methyladamantane.log \
         1_3_dimethyladamantane.molden 1_3_dimethyladamantane.log \
         1_3_5_trimethyladamantane.molden 1_3_5_trimethyladamantane.log; do
    if [ ! -f "$f" ]; then
        missing=$((missing + 1))
    fi
done

if [ "$missing" -gt 0 ]; then
    echo
    echo "ERROR: $missing file(s) missing."
    echo "Place them manually from https://github.com/AkimovLab/Project_adamantanes"
    exit 1
fi

echo
echo "All inputs present. Run 'bash run_all.sh'."
