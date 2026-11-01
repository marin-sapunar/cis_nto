#!/usr/bin/env python3
"""
Extract the leading single-excitation weight for state 1 from a CP2K
TDDFPT log.

The JPCL paper (doi:10.1021/acs.jpclett.5c01417) reports the "HOMO ->
LUMO weight" as the squared amplitude of the dominant single excitation
in the CIS expansion of the target state.  That is:

    w_leading = max_ia ( X_{i,a}^2 )

where X_{i,a} are the CP2K "Excitation amplitude" values printed in the
Excitation analysis block.  This is the quantity our NTO occupations
should be compared against, not the raw NTO eigenvalues (which are
related but not equal for multi-configurational states).

Usage:
    python3 extract_leading_weight.py <log_path> [--state N]

Where <log_path> is a CP2K out.log containing an "Excitation analysis"
block, and --state N selects which state to analyse (default: 1).
"""
import argparse
import re
import sys
from pathlib import Path


def parse_state_block(text, state):
    """Return list of (occ, vir, amplitude) for the given state.

    Returns None if the state header is not found.
    """
    # Locate the excitation analysis block
    m = re.search(r"Excitation analysis", text)
    if not m:
        return None
    body = text[m.end():]

    # Find the header line for our target state, e.g. "  1   2.27982 eV"
    pat = re.compile(r"^\s+" + str(state) + r"\s+[\d.]+\s+eV\s*$", re.MULTILINE)
    m = pat.search(body)
    if not m:
        return None

    after = body[m.end():]
    excitations = []
    for line in after.split("\n"):
        # Stop at the next state header
        if re.match(r"^\s+" + str(state + 1) + r"\s+[\d.]+\s+eV", line):
            break
        # Stop at the end of the block (blank line or a new section)
        if line.strip() == "":
            if excitations:
                break
            continue
        if not re.match(r"^\s+\d+\s+\d+\s+-?[\d.]+\s*$", line):
            if excitations:
                break
            continue
        parts = line.split()
        try:
            occ = int(parts[0])
            vir = int(parts[1])
            amp = float(parts[2])
        except (ValueError, IndexError):
            continue
        excitations.append((occ, vir, amp))
    return excitations


def leading_weight(log_path, state=1):
    """Return dict with leading weight, occ, vir, and total weight."""
    text = Path(log_path).read_text(errors="replace")
    ex = parse_state_block(text, state)
    if ex is None:
        return None
    if not ex:
        return None
    total_w = sum(a * a for _, _, a in ex)
    w, occ, vir = max(((a * a, o, v) for o, v, a in ex), key=lambda t: t[0])
    return {
        "leading_weight": w,
        "leading_occ": occ,
        "leading_vir": vir,
        "total_weight": total_w,
        "n_excitations": len(ex),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log_path", help="CP2K out.log with Excitation analysis block")
    ap.add_argument("--state", type=int, default=1,
                    help="state index (default 1)")
    args = ap.parse_args()

    r = leading_weight(args.log_path, state=args.state)
    if r is None:
        print(f"Could not find state {args.state} in {args.log_path}")
        sys.exit(1)

    label = Path(args.log_path).parent.name
    print(f"{label}: state {args.state}")
    print(f"  leading weight = {r['leading_weight']:.6f} "
          f"({r['leading_occ']} -> {r['leading_vir']})")
    print(f"  total weight   = {r['total_weight']:.6f}")
    print(f"  n_excitations  = {r['n_excitations']}")


if __name__ == "__main__":
    main()
