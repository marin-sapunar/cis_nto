"""Subprocess wrapper for cis_overlap.exe."""
import os
import subprocess
from pathlib import Path


def find_binary():
    candidates = [
        os.environ.get("CIS_OVERLAP_BIN"),
        "cis_overlap.exe",
        "./build/bin/cis_overlap.exe",
        "../build/bin/cis_overlap.exe",
        str(Path.home() / "cis_nto" / "build" / "bin" / "cis_overlap.exe"),
    ]
    for c in candidates:
        if c and Path(c).is_file():
            return str(c)
    raise FileNotFoundError("cis_overlap.exe not found. Set CIS_OVERLAP_BIN.")


def run_overlap_calculator(wf_1, wf_2, format="cp2k", nstates=3,
                           method="reprojection", binary=None):
    binary = binary or find_binary()
    cmd = [binary, "-in", format, "-out", "molden_cart", f"{wf_1}:{wf_2}"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"cis_overlap.exe failed:\n{result.stderr}")
    return result.stdout, result.stderr
