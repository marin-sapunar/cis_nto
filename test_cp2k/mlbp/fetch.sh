#!/usr/bin/env bash
#
# Fetch the ML-BP B3LYP test input from the AkimovLab repository.
#
# Source: https://github.com/AkimovLab/Project_MLBP_excitonic
# The upstream repository stores CP2K output in .tar.gz archives.
#
set -eu

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

BASE="https://raw.githubusercontent.com/AkimovLab/Project_MLBP_excitonic/main/step3"
TARBALL="B3LYP_50_states.tar.gz"
EXTRACTED_DIR="B3LYP_50states"
MOLDEN="b3lyp_opted-MOS-1_0.molden"
LOG="out.log"

if [ -f "$MOLDEN" ] && [ -f "$LOG" ]; then
    echo "Input already present:"
    ls -la "$MOLDEN" "$LOG"
    exit 0
fi

echo "Downloading $TARBALL ..."
curl -fsSL -o "$TARBALL" "$BASE/$TARBALL"

echo "Extracting ..."
tar -xzf "$TARBALL"

# Locate the two files
molden_src=$(find . -name "$MOLDEN" -print -quit)
log_src=$(find . -name "$LOG" -print -quit)

if [ -z "$molden_src" ] || [ -z "$log_src" ]; then
    echo "ERROR: expected files not found in tarball"
    echo "  Wanted: $MOLDEN and $LOG"
    exit 1
fi

mv "$molden_src" "$MOLDEN"
mv "$log_src" "$LOG"

# Clean up
rm -rf "$EXTRACTED_DIR"
rm -f "$TARBALL"

echo
echo "Done:"
ls -la "$MOLDEN" "$LOG"
echo
echo "Next: bash run.sh"
