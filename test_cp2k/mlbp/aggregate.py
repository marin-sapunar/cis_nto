#!/usr/bin/env python3
"""
Aggregate reconstruction, leading CI weight, and NTO occupation data
from all ML-BP functional runs into a single results.json.

Reads:
  nto_<FUNC>.log       - reader output (reconstruction notice + NTO
                         construction; the log also carries the CP2K
                         excitation amplitudes we extract here)
  nto_<FUNC>.molden    - NTO Molden (eigenvalue occupations)
  <FUNC>_50states/out.log - CP2K TDDFPT log with the "Excitation
                         analysis" block (source of the leading CI
                         amplitude)

Writes:
  results.json         - one record per functional

Comparison to the paper
-----------------------
The JPCL paper (doi:10.1021/acs.jpclett.5c01417, Table 1) reports
the **leading single-excitation weight** (max of X_ia^2 over all
excitations in the target state).  This is a property of the CI
vector, not of the NTO eigenvalue spectrum.  For a single-
configurational state the two coincide; for a multi-configurational
state they differ, sometimes substantially (e.g. BHandHLYP S1:
leading CI weight 0.8365 vs. leading NTO eigenvalue 0.9149).

Earlier versions of this script compared the NTO eigenvalue
(Sym=1001a) against the paper's leading CI weight, which produced
apparent deviations of 13-16 pp for range-separated and global
hybrids.  Those deviations were a quantity-definition mismatch, not
a numerical error: the reader had always been parsing the log
correctly.  The `delta_paper` field below is now computed against
the leading CI weight.
"""
import json
import re
from pathlib import Path
from datetime import datetime, timezone

HERE = Path(__file__).parent

# Paper reference values (Table 1 of J. Phys. Chem. Lett. 2024,
# doi:10.1021/acs.jpclett.5c01417), with state index.
#
# NOTE: for wB97X and wB97XD the paper's state is 38, not state 1.
# The state-1 leading weight is still reported (informational), but
# delta_paper remains null for those two functionals because the
# quantities are not directly comparable.
PAPER_REFS = {
    "PBEh0125":  {"state": 1,  "value": 0.993, "note": "HOMO -> LUMO"},
    "B3LYP":     {"state": 1,  "value": 0.979, "note": "HOMO -> LUMO"},
    "PBE0":      {"state": 1,  "value": 0.967, "note": "HOMO -> LUMO"},
    "BHandHLYP": {"state": 1,  "value": 0.837, "note": "HOMO -> LUMO"},
    "CAM_B3LYP": {"state": 1,  "value": 0.769, "note": "HOMO -> LUMO"},
    "wB97XD":    {"state": 38, "value": 0.185, "note": "HOMO-3 -> LUMO+5"},
    "wB97X":     {"state": 38, "value": 0.153, "note": "HOMO-5 -> LUMO+7"},
}

# Mapping from the FUNCTIONALS keys above to the on-disk subdirectory
# holding the CP2K log for that functional.
LOG_DIRS = {
    "PBEh0125":  "PBEh0125_50states",
    "B3LYP":     ".",
    "PBE0":      "PBE0_50_states",
    "BHandHLYP": "BHandHLYP_50states",
    "CAM_B3LYP": "CAM_B3LYP_50_states",
    "wB97XD":    "wb97XD_50_staes",
    "wB97X":     "wB97X_50_states",
}

FUNCTIONALS = list(PAPER_REFS.keys())


def parse_reconstruction(log_path):
    """Extract the reconstruction notice from nto_<FUNC>.log."""
    if not log_path.exists():
        return None
    text = log_path.read_text(errors="replace")
    rec = {}
    m = re.search(r"log references\s*:\s*(\d+)", text)
    if m: rec["log_references"] = int(m.group(1))
    m = re.search(r"molden contains\s*:\s*(\d+)", text)
    if m: rec["molden_contains"] = int(m.group(1))
    m = re.search(r"reconstruction method\s*:\s*(\S+)", text)
    if m: rec["method"] = m.group(1)
    m = re.search(r"dropping\s+(\d+)\s+virtuals", text)
    if m: rec["n_dropped"] = int(m.group(1))
    return rec or None


def parse_dominant_nto(molden_path):
    """Return NTO eigenvalues from the Molden file.

    The first two NTOs written by cis_nto are conventionally labelled
    1001a and 1002a; their occupations are the leading singular values
    of the SVD of the CI amplitude matrix.  These are NOT the same as
    the leading single-excitation weight (see module docstring).
    """
    if not molden_path.exists():
        return {}
    occ_by_sym = {}
    current_sym = None
    with molden_path.open(errors="replace") as f:
        for line in f:
            m = re.match(r"^\s*Sym=\s*(\S+)", line)
            if m:
                current_sym = m.group(1)
                continue
            m = re.match(r"^\s*Occup=\s*([\d.eE+-]+)", line)
            if m and current_sym is not None:
                try:
                    occ_by_sym[current_sym] = float(m.group(1))
                except ValueError:
                    pass
                current_sym = None
    return {
        "leading_nto_eigenvalue": occ_by_sym.get("1001a"),
        "second_nto_eigenvalue":  occ_by_sym.get("1002a"),
    }


def parse_leading_ci_weight(cp2k_log_path, state=1):
    """Return leading single-excitation weight max_ia (X_ia^2).

    Reads the CP2K TDDFPT "Excitation analysis" block.  This is the
    quantity reported by the JPCL paper Table 1.
    """
    if not cp2k_log_path.exists():
        return None
    text = cp2k_log_path.read_text(errors="replace")

    m = re.search(r"Excitation analysis", text)
    if not m:
        return None
    body = text[m.end():]

    pat = re.compile(r"^\s+" + str(state) + r"\s+[\d.]+\s+eV\s*$", re.MULTILINE)
    m = pat.search(body)
    if not m:
        return None

    after = body[m.end():]
    best = (0.0, None, None)   # (amp^2, occ, vir)
    total = 0.0
    n_exc = 0
    for line in after.split("\n"):
        if re.match(r"^\s+" + str(state + 1) + r"\s+[\d.]+\s+eV", line):
            break
        if line.strip() == "":
            if n_exc:
                break
            continue
        if not re.match(r"^\s+\d+\s+\d+\s+-?[\d.]+\s*$", line):
            if n_exc:
                break
            continue
        parts = line.split()
        try:
            occ = int(parts[0])
            vir = int(parts[1])
            amp = float(parts[2])
        except (ValueError, IndexError):
            continue
        total += amp * amp
        n_exc += 1
        if amp * amp > best[0]:
            best = (amp * amp, occ, vir)

    if n_exc == 0:
        return None
    return {
        "leading_weight": best[0],
        "leading_occ": best[1],
        "leading_vir": best[2],
        "total_weight": total,
        "n_excitations": n_exc,
    }


def main():
    results = {
        "schema_version": "1.1",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dataset": "ML-BP (AkimovLab/Project_MLBP_excitonic, step3)",
        "paper_doi": "10.1021/acs.jpclett.5c01417",
        "notes": (
            "delta_paper is computed against the leading single-excitation "
            "weight (max_ia X_ia^2), the quantity reported by the paper. "
            "NTO eigenvalues are reported as supplemental information and "
            "are not directly comparable to the paper."
        ),
        "functionals": [],
    }

    for func in FUNCTIONALS:
        log = HERE / f"nto_{func}.log"
        molden = HERE / f"nto_{func}.molden"
        rec = parse_reconstruction(log)
        nto = parse_dominant_nto(molden)
        paper = PAPER_REFS[func]

        # Locate the CP2K log for the leading-CI-weight extraction
        log_dir = LOG_DIRS.get(func)
        cp2k_log = HERE / log_dir / "out.log" if log_dir else None
        lead = parse_leading_ci_weight(cp2k_log, state=1) if cp2k_log else None

        record = {
            "functional": func,
            "state": 1,
            "reconstruction": rec,
            "leading_ci_weight": lead,
            "nto": nto,
            "paper": paper,
            "delta_paper": None,
        }

        # Compute the paper deviation only when the paper's state matches
        # the state we extracted (both must be 1).
        if (lead is not None
                and paper["state"] == 1
                and lead.get("leading_weight") is not None):
            record["delta_paper"] = round(
                lead["leading_weight"] - paper["value"], 6
            )

        results["functionals"].append(record)

    out = HERE / "results.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"Wrote {out}")
    print(f"  {len(results['functionals'])} functionals")

    # Summary table
    print()
    hdr = f"{'functional':<12} {'state':>5}  {'lead W':>8}  {'paper':>8}  {'delta':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in results["functionals"]:
        lead = r["leading_ci_weight"]
        lead_w = lead["leading_weight"] if lead else None
        state = r["state"]
        paper = r["paper"]
        delta = r["delta_paper"]
        lead_tag = f"{lead_w:.4f}" if lead_w is not None else "  n/a "
        paper_tag = f"{paper['value']:.4f}" if paper["state"] == state else " (s38)"
        delta_tag = f"{delta:+.4f}" if delta is not None else "  n/a "
        print(f"{r['functional']:<12} {state:>5}  {lead_tag:>8}  {paper_tag:>8}  {delta_tag:>8}")


if __name__ == "__main__":
    main()
