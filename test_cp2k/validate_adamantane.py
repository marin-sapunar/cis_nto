#!/usr/bin/env python3
"""
Adamantane family capability validator for cis_nto.

Demonstrates the CP2K reader on a chemically meaningful molecular
series: adamantane, 1-methyl-, 1,3-dimethyl-, 1,3,5-trimethyl-
adamantane.

Usage:
    python3 validate_adamantane.py [--exe PATH]
Exit code 0 if all checks pass.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (input_name, output_name, formula)
# Inputs use underscores (e.g. 1_methyladamantane.molden)
# Outputs strip them (e.g. s1_virt_1methyladamantane.cube)
FAMILY = [
    ("adamantane",                 "adamantane",                "C10H16"),
    ("1_methyladamantane",         "1methyladamantane",         "C11H18"),
    ("1_3_dimethyladamantane",     "13dimethyladamantane",      "C12H20"),
    ("1_3_5_trimethyladamantane",  "135trimethyladamantane",    "C13H22"),
]

# Helper: input pair for reader (molden, log)
def input_paths(prefix):
    return HERE / f"{prefix}.molden", HERE / f"{prefix}.log"

# Helper: cube path for a given kind and output name
def cube_path(kind, out_name):
    return HERE / "cubes" / f"{kind}_{out_name}.cube"

# Helper: NTO molden path
def nto_molden_path(out_name):
    return HERE / f"nto_{out_name}.molden"

# Helper: NTO log path
def nto_log_path(out_name):
    return HERE / f"nto_{out_name}.log"



def banner(t):
    print(); print("=" * 68); print(t); print("=" * 68)


def report(label, passed, detail=""):
    mark = "PASS" if passed else "FAIL"
    line = f"  [{mark}] {label}"
    if detail:
        line += f"  --  {detail}"
    print(line)
    return passed


def check_inputs():
    banner("Section 1 -- Input availability")
    print("  Adamantane inputs are fetched at run time by fetch.sh.")
    print()
    ok = True
    for in_name, out_name, formula in FAMILY:
        m, l = input_paths(in_name)
        have = m.exists() and l.exists()
        report(f"{in_name:<28} inputs present", have, f"formula {formula}")
        ok &= have
    if not ok:
        print()
        print("  --> run:  bash fetch.sh")
    return ok


def run_reader(exe, in_name):
    molden, log = input_paths(in_name)
    if not (molden.exists() and log.exists()):
        return False, "", "inputs missing"
    try:
        r = subprocess.run(
            [str(exe), "-in", "cp2k", f"{molden}:{log}"],
            capture_output=True, timeout=180, text=True,
            cwd=str(HERE),
        )
        return r.returncode == 0, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return False, "", "timeout"


def check_reader(exe):
    banner("Section 2 -- Reader invocation on the full family")
    ok = True
    for in_name, out_name, _ in FAMILY:
        good, out, err = run_reader(exe, in_name)
        combined = (out + err).lower()
        clean = good and "error" not in combined
        report(f"{in_name:<28} reader runs", clean, f"exit={0 if good else 1}")
        ok &= clean
    return ok


def check_cubes():
    banner("Section 3 -- NTO cube generation")
    cube_dir = HERE / "cubes"
    if not cube_dir.is_dir():
        report("cubes/ directory present", False, "run run_all.sh first")
        return False
    ok = True
    print(f"  {'Molecule':<28} {'homo':<6} {'lumo':<6} {'s1_occ':<8} {'s1_virt':<8}")
    print("  " + "-" * 60)
    for in_name, out_name, _ in FAMILY:
        have = {}
        for kind in ("homo", "lumo", "s1_occ", "s1_virt"):
            p = cube_path(kind, out_name)
            have[kind] = p.exists()
        all_ok = all(have.values())
        cells = " ".join("Y" if have[k] else "n" for k in
                         ("homo", "lumo", "s1_occ", "s1_virt"))
        report(f"{in_name:<28} cubes", all_ok, cells)
        ok &= all_ok
    return ok


def parse_cube_stats(cube_path):
    """Return (n_nonzero, l2_integral) for a Gaussian cube file."""
    with cube_path.open(errors="replace") as f:
        f.readline(); f.readline()
        parts = f.readline().split()
        natoms = int(parts[0])
        voxels = []
        for _ in range(3):
            p = f.readline().split()
            voxels.append((int(p[0]),
                           (float(p[1]), float(p[2]), float(p[3]))))
        for _ in range(natoms):
            f.readline()
        n_nonzero = 0
        l2_sq = 0.0
        for line in f:
            for p in line.split():
                try:
                    v = float(p)
                except ValueError:
                    continue
                if abs(v) > 1e-8:
                    n_nonzero += 1
                l2_sq += v * v
    d1 = voxels[0][1]; d2 = voxels[1][1]; d3 = voxels[2][1]
    cross = (d2[1]*d3[2] - d2[2]*d3[1],
             d2[2]*d3[0] - d2[0]*d3[2],
             d2[0]*d3[1] - d2[1]*d3[0])
    dV = abs(d1[0]*cross[0] + d1[1]*cross[1] + d1[2]*cross[2])
    return n_nonzero, (l2_sq ** 0.5) * (dV ** 0.5)


def check_physics():
    banner("Section 4 -- Physics consistency checks on NTO cubes")
    print("  L2 integral ~ 1.0; nonzero-cell count above a floor.")
    print()
    ok = True
    print(f"  {'Molecule':<28} {'L2':<8} {'nonzero':<10} {'verdict':<8}")
    print("  " + "-" * 60)
    for in_name, out_name, _ in FAMILY:
        cube = cube_path("s1_virt", out_name)
        if not cube.exists():
            continue
        n_nonzero, l2 = parse_cube_stats(cube)
        l2_ok      = abs(l2 - 1.0) < 5e-2
        nonzero_ok = n_nonzero > 100
        verdict = "OK" if (l2_ok and nonzero_ok) else "CHECK"
        print(f"  {in_name:<28} {l2:<8.4f} {n_nonzero:<10d} {verdict:<8}")
        ok &= (l2_ok and nonzero_ok)
    return ok


def check_methylation_trend():
    banner("Section 5 -- Chemical trend: methylation changes NTO character")
    values = {}
    for in_name, out_name, _ in FAMILY:
        log = nto_molden_path(out_name)  # Occup= lives in the molden
        if not log.exists():
            continue
        body = log.read_text(errors="replace")
        m = re.search(r"Occup=\s*([\d.eE+-]+)", body)
        if m:
            values[in_name] = float(m.group(1))

    if len(values) < 2:
        print("  (insufficient data; run run_all.sh first)")
        return False

    print(f"  {'Molecule':<28} {'S1 dominant NTO':<16}")
    print("  " + "-" * 46)
    for in_name, out_name, _ in FAMILY:
        if in_name in values:
            print(f"  {in_name:<28} {values[in_name]:<16.6f}")

    first = list(values.values())[0]
    last  = list(values.values())[-1]
    change = abs(last - first)
    sensible = change > 1e-3
    print()
    report(f"endpoint change > 1e-3 ({change:.4f})", sensible)
    return sensible


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", type=Path,
                    default=HERE.parent / "build" / "bin" / "cis_nto.exe")
    args = ap.parse_args()
    print()
    print("Adamantane family capability validation")
    print(f"  binary: {args.exe}")
    print(f"  data:   {HERE}")
    results = [
        ("1. Input availability",  check_inputs()),
        ("2. Reader invocation",   check_reader(args.exe)),
        ("3. NTO cube generation", check_cubes()),
        ("4. Physics checks",      check_physics()),
        ("5. Methylation trend",   check_methylation_trend()),
    ]
    banner("Summary")
    all_pass = True
    for name, ok in results:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        all_pass &= ok
    print()
    print("  ALL ADAMANTANE CHECKS PASSED" if all_pass
          else "  SOME CHECKS FAILED")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
