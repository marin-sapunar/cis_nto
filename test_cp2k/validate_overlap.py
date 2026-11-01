#!/usr/bin/env python3
"""Self-contained overlap pipeline validation — produces T1/T2/T3/T5/T6 tables."""
import os, sys, subprocess, time
from pathlib import Path
import numpy as np

VAC = Path(os.environ.get('VAC', Path.home() / 'Project_ML_BP_revision1/03-OVERLAPS/di_vacancy'))
CIS_OVERLAP = Path(os.environ.get('CIS_BIN', Path.home() / 'cis_nto/build/bin')) / 'cis_overlap.exe'
OUT = VAC / f'py_validation_{time.strftime("%Y%m%d_%H%M%S")}'
OUT.mkdir(exist_ok=True, parents=True)

print(f"Binary: {CIS_OVERLAP}")
print(f"Data:   {VAC}")
print(f"Output: {OUT}\n")

if not CIS_OVERLAP.exists():
    sys.exit(f"FATAL: {CIS_OVERLAP} not found")

def load_swf(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        M = np.zeros((n, n))
        for i in range(n):
            M[i] = [float(x) for x in f.readline().split()]
    return M

def normalize_with_sign(M):
    d = np.diag(M)
    sign = np.where(d >= 0, 1.0, -1.0)
    d_abs = np.abs(d)
    d_abs = np.where(d_abs < 1e-15, 1.0, d_abs)
    return M / np.sqrt(np.outer(d_abs, d_abs)) * np.outer(sign, sign)

def run_overlap(m1, l1, m2, l2, tnex=1, timeout=90):
    swf = VAC / 's_wf'
    if swf.exists(): swf.unlink()
    cmd = [str(CIS_OVERLAP), '-in', 'cp2k', '-alg', 'NTO', '-tnex', str(tnex),
           f'{m1}:{l1}', f'{m2}:{l2}']
    try:
        r = subprocess.run(cmd, cwd=VAC, timeout=timeout,
                           capture_output=True, text=True)
        if not swf.exists():
            return None, f"no s_wf (rc={r.returncode}): {r.stderr[:200]}"
        return load_swf(swf), None
    except subprocess.TimeoutExpired:
        return None, 'timeout'

# --- Discover complete frames ---
frames = []
for n in range(1, 100):
    m = VAC / f'p4-di-PBE_cont_{n}-MOS-1_0.molden'
    l = VAC / f'run_cont_{n}.log'
    if not (m.exists() and l.exists()): break
    try:
        with open(l) as f:
            if 'TDDFPT run converged' not in f.read(): break
    except: break
    frames.append(n)

if len(frames) < 2:
    sys.exit(f"FATAL: need >=2 frames, found {frames}")

print(f"Frames: {frames}  total={len(frames)}\n")
BAR = "=" * 78
THIN = "-" * 78

# ============================================================
# T1: Self-overlap identity (-tnex 3)
# ============================================================
print(f"{BAR}\n T1 — SELF-OVERLAP IDENTITY (-tnex 3)\n{BAR}")
print(f"{'frame':>5} | {'GS':>7} {'S1':>7} {'S2':>7} {'S3':>7} {'S4':>7} {'S5':>7} | status")
print(THIN)
t1_pass = 0
for n in frames:
    m = f'p4-di-PBE_cont_{n}-MOS-1_0.molden'
    l = f'run_cont_{n}.log'
    M, err = run_overlap(m, l, m, l, tnex=3)
    if M is None:
        print(f"{n:>5} | ERROR: {err}")
        continue
    Mn = normalize_with_sign(M)  # <--- Added normalization here
    np.save(OUT / f't1_self_{n}.npy', Mn)
    d = np.abs(np.diag(Mn))[:6]
    ok = abs(d[0] - 1) < 1e-3 and all(abs(e - 1) < 1e-3 for e in d[1:])
    status = "PASS" if ok else "FAIL"
    if ok: t1_pass += 1
    print(f"{n:>5} | " + " ".join(f"{x:>7.5f}" for x in d) + f" | {status}")
print(THIN)
print(f"T1: {t1_pass}/{len(frames)} PASS\n")

# ============================================================
# T2: Consecutive tracking (-tnex 1, normalized)
# ============================================================
print(f"{BAR}\n T2 — CONSECUTIVE TRACKING (S1 row off-diagonals)\n{BAR}")
print(f"{'pair':>8} | {'S1↔S2':>8} {'S1↔S3':>8} {'S1↔S4':>8} {'S1↔S5':>8}")
print(THIN)
t2_pairs = []
for k in range(len(frames) - 1):
    a, b = frames[k], frames[k+1]
    M, err = run_overlap(f'p4-di-PBE_cont_{a}-MOS-1_0.molden', f'run_cont_{a}.log',
                         f'p4-di-PBE_cont_{b}-MOS-1_0.molden', f'run_cont_{b}.log', tnex=1)
    if M is None:
        print(f"{a}→{b:>4} | ERROR: {err}")
        continue
    Mn = normalize_with_sign(M)
    np.save(OUT / f't2_cross_{a}_{b}.npy', Mn)
    row = np.abs(Mn[1, 2:6])
    t2_pairs.append((a, b, row))
    print(f"{a}→{b:<4} | " + " ".join(f"{x:>8.5f}" for x in row))
print(THIN)
print(f"T2: {len(t2_pairs)}/{len(frames)-1} pairs\n")

# ============================================================
# T3: Reference S(1→N) (-tnex 1, normalized)
# ============================================================
print(f"{BAR}\n T3 — REFERENCE DECAY S(1→N) (S1 row off-diagonals)\n{BAR}")
first = frames[0]
print(f"{'pair':>8} | {'S1↔S2':>8} {'S1↔S3':>8} {'S1↔S4':>8} {'S1↔S5':>8}")
print(THIN)
t3_pairs = []
for k in range(1, len(frames)):
    b = frames[k]
    M, err = run_overlap(f'p4-di-PBE_cont_{first}-MOS-1_0.molden', f'run_cont_{first}.log',
                         f'p4-di-PBE_cont_{b}-MOS-1_0.molden', f'run_cont_{b}.log', tnex=1)
    if M is None:
        print(f"{first}→{b:>4} | ERROR: {err}")
        continue
    Mn = normalize_with_sign(M)
    np.save(OUT / f't3_ref_{first}_{b}.npy', Mn)
    row = np.abs(Mn[1, 2:6])
    t3_pairs.append((b - first, row))
    print(f"{first}→{b:<4} | " + " ".join(f"{x:>8.5f}" for x in row))
print(THIN)
print(f"T3: {len(t3_pairs)} points\n")

# ============================================================
# T5 + T6: rate fit + classification
# ============================================================
if len(t3_pairs) >= 3:
    dts = np.array([p[0] for p in t3_pairs], dtype=float)
    cols = {
        'S1↔S2': np.array([p[1][0] for p in t3_pairs]),
        'S1↔S3': np.array([p[1][1] for p in t3_pairs]),
        'S1↔S4': np.array([p[1][2] for p in t3_pairs]),
        'S1↔S5': np.array([p[1][3] for p in t3_pairs]),
    }

    print(f"{BAR}\n T5 — DECAY RATE FIT ({len(dts)} points, Δt={int(dts[0])}..{int(dts[-1])})\n{BAR}")
    print(f"{'column':<10} | {'slope (×10⁻³/Δt)':>18} | {'R²':>6} | {'mean':>8} | {'σ':>8} | {'CV':>6}")
    print(THIN)
    fit_results = {}
    for name, arr in cols.items():
        slope, intercept = np.polyfit(dts, arr, 1)
        pred = slope * dts + intercept
        ss_res = np.sum((arr - pred)**2)
        ss_tot = np.sum((arr - arr.mean())**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
        cv = arr.std() / arr.mean() if arr.mean() > 0 else 0
        fit_results[name] = dict(slope=slope, r2=r2, cv=cv, mean=arr.mean(), std=arr.std())
        print(f"{name:<10} | {slope*1e3:>18.3f} | {r2:>6.3f} | {arr.mean():>8.5f} | {arr.std():>8.5f} | {cv:>6.3f}")
    print(THIN + "\n")

    print(f"{BAR}\n T6 — PHYSICS vs RESIDUAL CLASSIFICATION\n{BAR}")
    print(f"  Rules: CV<0.05 → RESIDUAL | CV>0.10 + monotone → PHYSICS | CV>0.20 → PHYSICS(noisy) | else MIXED")
    print(THIN)
    print(f"{'column':<10} | {'CV':>6} | {'monotone':>9} | {'verdict':<14} | interpretation")
    print(THIN)
    phys_cols, resid_cols = [], []
    for name, arr in cols.items():
        info = fit_results[name]
        cv = info['cv']
        diffs = np.diff(arr)
        mono_up = all(d > -1e-3 for d in diffs)
        mono_down = all(d < 1e-3 for d in diffs)
        if cv < 0.05:
            verdict = "RESIDUAL"; resid_cols.append(name)
            interp = f"constant ≈ CP2K [5D7F] residual"
        elif cv > 0.20:
            verdict = "PHYSICS(noisy)"; phys_cols.append(name)
            interp = f"fluctuating mixing, slope {info['slope']*1e3:+.2f}e-3"
        elif cv > 0.10 and (mono_up or mono_down):
            verdict = "PHYSICS"; phys_cols.append(name)
            interp = f"state mixing, slope {info['slope']*1e3:+.2f}e-3"
        else:
            verdict = "MIXED"
            interp = "weak trend / noise"
        mono = "↑" if mono_up else "↓" if mono_down else "no"
        print(f"{name:<10} | {cv:>6.3f} | {mono:>9} | {verdict:<14} | {interp}")
    print(THIN)
    print(f"\n  PHYSICS:  {', '.join(phys_cols) if phys_cols else '(none)'}")
    print(f"  RESIDUAL: {', '.join(resid_cols) if resid_cols else '(none)'}")
    print(f"\n  Note: at {len(dts)} points, CV uncertainty ≈ ±{20/np.sqrt(len(dts)):.0f}%")
    print(f"  Rerun with ≥15 points for tighter classification.\n")

print(BAR)
print(f" DONE — artifacts in {OUT}/")
print(BAR)
