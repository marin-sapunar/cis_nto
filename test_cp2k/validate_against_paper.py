#!/usr/bin/env python3
"""
Parse a CP2K TDDFPT log's 'Excitation analysis' block and compare the
leading HOMO->LUMO transition amplitude for a chosen state against a
user-supplied reference value.

Usage:
    validate_against_paper.py <cp2k_log> <state> [ref_weight]

<cp2k_log> is the CP2K output file (NOT the cis_nto banner; those are
written as nto_<tag>.log and do not contain the excitation analysis).
<state> is 1-based. If [ref_weight] is given, the script checks the
leading HOMO->LUMO amplitude against it to 1e-3 and exits nonzero on
mismatch. If omitted, the script just prints the parsed transitions.

Reference values must be supplied by the user. They depend on the
system and on the quantity the user wants to compare — the CP2K log
reports raw CI coefficients, which are not necessarily identical to
the state-averaged weights printed in a paper's Table 1.

Exit code: 0 on success, 1 on mismatch, 2 on parse failure.
"""
import re
import sys


def parse_excitations(logfile):
    """Return {state: [(occ, virt, amp), ...]} from the excitation analysis."""
    states = {}
    current = None
    in_block = False
    with open(logfile) as fh:
        for line in fh:
            s = line.rstrip()
            if "Excitation analysis" in s:
                in_block = True
                continue
            if not in_block:
                continue
            # End-of-block heuristics
            if "TDDFPT" in s and "analysis" not in s:
                in_block = False
                continue
            m = re.match(r"\s*(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s*$", s)
            if m:
                occ, virt, amp = int(m.group(1)), int(m.group(2)), float(m.group(3))
                if current is not None:
                    states.setdefault(current, []).append((occ, virt, amp))
                continue
            m = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)", s)
            if m:
                current = int(m.group(1))
                states.setdefault(current, [])
    return states


def leading_transition(transitions, occ_ref=160, virt_ref=161):
    """Return the amplitude of the dominant HOMO->LUMO (occ_ref, virt_ref)."""
    best = 0.0
    for occ, virt, amp in transitions:
        if occ == occ_ref and virt == virt_ref:
            return abs(amp)
        best = max(best, abs(amp))
    return best


def main():
    if len(sys.argv) < 3:
        print("Usage: validate_against_paper.py <log> <state> [weight] [f_osc]")
        sys.exit(2)
    logfile = sys.argv[1]
    state = int(sys.argv[2])
    ref_weight = float(sys.argv[3]) if len(sys.argv) > 3 else None
    ref_fosc   = float(sys.argv[4]) if len(sys.argv) > 4 else None

    states = parse_excitations(logfile)
    if not states:
        print(f"FAIL: no excitation analysis found in {logfile}")
        sys.exit(2)

    if state not in states:
        print(f"FAIL: state {state} not found (have {sorted(states)[:5]}...)")
        sys.exit(2)

    trans = states[state]
    print(f"State {state}: {len(trans)} transitions")

    # Leading contribution
    lead = max(trans, key=lambda t: abs(t[2]))
    print(f"  leading: occ={lead[0]} virt={lead[1]} amp={lead[2]:.6f}")

    # HOMO->LUMO specifically (occ=160, virt=161)
    hl_amp = leading_transition(trans)
    print(f"  HOMO->LUMO amplitude: {hl_amp:.6f}")

    ok = True
    if ref_weight is not None:
        d = abs(hl_amp - ref_weight)
        if d <= 1e-3:
            print(f"  PASS  weight {hl_amp:.4f} vs ref {ref_weight:.4f}")
        else:
            print(f"  FAIL  weight {hl_amp:.4f} vs ref {ref_weight:.4f} (d={d:.4f})")
            ok = False

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
