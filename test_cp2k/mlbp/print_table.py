import re
from pathlib import Path

HERE = Path(__file__).parent
PAPER = {
    "PBEh0125":  (1, 0.993), "B3LYP": (1, 0.979), "PBE0": (1, 0.967),
    "BHandHLYP": (1, 0.837), "CAM_B3LYP": (1, 0.769),
    "wB97XD":    (38, 0.185), "wB97X": (38, 0.153),
}
LOG_DIRS = {
    "PBEh0125":  "PBEh0125_50states", "B3LYP": "B3LYP_50_states",
    "PBE0":      "PBE0_50_states",    "BHandHLYP": "BHandHLYP_50states",
    "CAM_B3LYP": "CAM_B3LYP_50_states", "wB97XD": "wb97XD_50_staes",
    "wB97X":     "wB97X_50_states",
}

print(f"{'functional':<12} {'state':>5}  {'lead W':>8}  {'paper':>8}  {'delta':>8}")
print("-" * 50)

for func, (state, paper_val) in PAPER.items():
    log_dir = LOG_DIRS.get(func)
    log_path = HERE / log_dir / "out.log" if log_dir else None
    
    lead_w = None
    if log_path and log_path.exists():
        text = log_path.read_text(errors="replace")
        m = re.search(r"Excitation analysis", text)
        if m:
            body = text[m.end():]
            pat = re.compile(r"^\s+" + str(state) + r"\s+[\d.]+\s+eV\s*$", re.MULTILINE)
            m = pat.search(body)
            if m:
                after = body[m.end():]
                best = 0.0
                for line in after.split("\n"):
                    if re.match(r"^\s+" + str(state + 1) + r"\s+[\d.]+\s+eV", line): break
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            amp = float(parts[2])
                            best = max(best, amp * amp)
                        except ValueError: pass
                lead_w = best

    if lead_w is not None:
        delta = lead_w - paper_val
        print(f"{func:<12} {state:>5}  {lead_w:>8.4f}  {paper_val:>8.4f}  {delta:>+8.4f}")
    else:
        print(f"{func:<12} {state:>5}  {'n/a':>8}  {paper_val:>8.4f}  {'n/a':>8}")
