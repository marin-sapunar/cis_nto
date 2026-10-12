#!/usr/bin/env python3
"""
Aggregate reconstruction and NTO data from all ML-BP functional runs
into a single results.json.

Reads:
  nto_<FUNC>.log       - reader output (reconstruction notice)
  nto_<FUNC>.molden    - NTO Molden (occupations)

Writes:
  results.json         - one record per functional
"""
import json
import re
from pathlib import Path
from datetime import datetime, timezone

HERE = Path(__file__).parent

# Paper reference values (Table 1 of J. Phys. Chem. Lett. 2024,
# doi:10.1021/acs.jpclett.5c01417), with state index.
PAPER_REFS = {
    "PBEh0125":  {"state": 1, "value": 0.993, "note": "HOMO -> LUMO"},
    "B3LYP":     {"state": 1, "value": 0.979, "note": "HOMO -> LUMO"},
    "PBE0":      {"state": 1, "value": 0.967, "note": "HOMO -> LUMO"},
    "BHandHLYP": {"state": 1, "value": 0.837, "note": "HOMO -> LUMO"},
    "CAM_B3LYP": {"state": 1, "value": 0.769, "note": "HOMO -> LUMO"},
    "wB97XD":    {"state": 38, "value": 0.185, "note": "HOMO-3 -> LUMO+5"},
    "wB97X":     {"state": 38, "value": 0.153, "note": "HOMO-5 -> LUMO+7"},
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
    """
    Return the occupation of the dominant occupied NTO (Sym=2001a) and
    the dominant virtual NTO (Sym=1001a).
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
        "dominant_occupied_nto": occ_by_sym.get("2001a"),
        "dominant_virtual_nto":  occ_by_sym.get("1001a"),
    }


def main():
    results = {
        "schema_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dataset": "ML-BP (AkimovLab/Project_MLBP_excitonic, step3)",
        "paper_doi": "10.1021/acs.jpclett.5c01417",
        "functionals": [],
    }

    for func in FUNCTIONALS:
        log = HERE / f"nto_{func}.log"
        molden = HERE / f"nto_{func}.molden"
        rec = parse_reconstruction(log)
        occ = parse_dominant_nto(molden)
        paper = PAPER_REFS[func]

        record = {
            "functional": func,
            "state": 1,
            "reconstruction": rec,
            "nto": occ,
            "paper": paper,
            "delta_paper": None,
        }
        if occ.get("dominant_occupied_nto") is not None and paper["state"] == 1:
            record["delta_paper"] = round(
                occ["dominant_occupied_nto"] - paper["value"], 6
            )
        results["functionals"].append(record)

    out = HERE / "results.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"Wrote {out}")
    print(f"  {len(results['functionals'])} functionals")
    for r in results["functionals"]:
        occ = r["nto"].get("dominant_occupied_nto")
        delta = r["delta_paper"]
        tag = f"{occ:.4f}" if occ is not None else "  n/a "
        dtag = f"{delta:+.4f}" if delta is not None else "  n/a "
        print(f"  {r['functional']:<11} occ={tag}  delta_vs_paper={dtag}")


if __name__ == "__main__":
    main()
