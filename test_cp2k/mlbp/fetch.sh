#!/usr/bin/env bash
#
# Fetch an ML-BP test input from the AkimovLab repository.
#
# Usage:
#   bash fetch.sh              # default: B3LYP
#   bash fetch.sh PBEh0125     # any of the seven upstream tarballs
#
# Source: https://github.com/AkimovLab/Project_MLBP_excitonic
#
set -eu

FUNC="${1:-B3LYP}"

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

BASE="https://raw.githubusercontent.com/AkimovLab/Project_MLBP_excitonic/main/step3"
TARBALL="${FUNC}_50_states.tar.gz"

# Mapping from tarball functional name to the directory CP2K's tarball
# actually creates. Some entries differ in case, use an underscore, or
# contain a typo in the upstream tarball.
case "$FUNC" in
    B3LYP)      DIRNAME="." ;;
    PBEh0125)   DIRNAME="PBEh0125_50states" ;;
    BHandHLYP)  DIRNAME="BHandHLYP_50states" ;;
    PBE0)       DIRNAME="PBE0_50_states" ;;
    CAM_B3LYP)  DIRNAME="CAM_B3LYP_50_states" ;;
    wB97X)      DIRNAME="wB97X_50_states" ;;
    wB97XD)     DIRNAME="wb97XD_50_staes" ;;
    *)          DIRNAME="${FUNC}_50_states" ;;
esac

# Check if the directory for this functional already exists
existing=""
if [ -d "$DIRNAME" ]; then
    existing="$DIRNAME/"
fi

if [ -n "$existing" ] && [ -f "${existing}out.log" ]; then
    echo "Already extracted: $existing"
    ls -la "$existing"
    exit 0
fi

echo "Downloading $TARBALL ..."
curl -fsSL -o "$TARBALL" "$BASE/$TARBALL"

echo "Extracting ..."
tar -xzf "$TARBALL"
rm -f "$TARBALL"

# Find the directory the tarball extracted to
extracted=""
for d in ${FUNC}*/ ; do
    if [ -d "$d" ] && [ -f "${d}out.log" ]; then
        extracted="$d"
        break
    fi
done

if [ -z "$extracted" ]; then
    echo "ERROR: could not find extracted directory for $FUNC"
    echo "Contents of current directory:"
    ls -la
    exit 1
fi

echo
echo "Done. Files in $extracted:"
ls -la "$extracted"
echo
echo "Next: run 'bash run.sh $FUNC' (after updating the case block)."
