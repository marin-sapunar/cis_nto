"""
Compute S_11(i -> i+1) for all adjacent-frame pairs.

Uses the CIS determinant-overlap formula:
  <Φ_i^a | Φ_j^b> = <i|j><a|b> - <i|b><a|j>
applied through the MO cross-overlap S_MO and the sparse CI coefficients.

For a closed-shell singlet CIS, the diagonal (state 1 to state 1) overlap is:
  S_11 = 2*T1 - T2
where
  T1 = sum_{i,j,a,b} X_i[i,a] S_oo[i,j] X_j[j,b] S_vv[a,b]
  T2 = sum_{i,j,a,b} X_i[i,a] S_ov[i,b] X_j[j,b] S_vo[a,j]

(Triplets use T1 - T2; singlets use 2*T1 - T2.)
"""
import re
import numpy as np

def load_swf(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        return np.loadtxt(f, max_rows=n)

def parse_ci_state(log_file, state_num=1):
    with open(log_file) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if "Excitation analysis" in l)
    state_re = re.compile(r'^\s+' + str(state_num) + r'\s+\d+\.\d+ eV')
    next_re = re.compile(r'^\s+\d+\s+\d+\.\d+ eV')
    i = start
    while i < len(lines) and not state_re.match(lines[i]):
        i += 1
    if i >= len(lines):
        return []
    i += 1
    ci = []
    while i < len(lines) and not next_re.match(lines[i]):
        p = lines[i].split()
        if len(p) == 3:
            try:
                ci.append((int(p[0]), int(p[1]), float(p[2])))
            except ValueError:
                pass
        i += 1
    return ci

def build_X(ci, n_occ=155, n_virt=50):
    X = np.zeros((n_occ, n_virt))
    for occ, virt, amp in ci:
        v = virt - (n_occ + 1)
        if 0 <= v < n_virt:
            X[occ-1, v] = amp
    return X

def state_overlap(X_i, X_j, S_MO, n_occ=155, singlet=True):
    """CIS state overlap (spin-adapted)."""
    S_oo = S_MO[:n_occ, :n_occ]
    S_ov = S_MO[:n_occ, n_occ:]
    S_vo = S_MO[n_occ:, :n_occ]
    S_vv = S_MO[n_occ:, n_occ:]

    # T1 = sum_{i,j,a,b} X_i[i,a] S_oo[i,j] X_j[j,b] S_vv[a,b]
    B = X_i.T @ S_oo @ X_j
    T1 = (B * S_vv).sum()
    # T2 = sum_{i,j,a,b} X_i[i,a] S_ov[i,b] X_j[j,b] S_vo[a,j]
    U = X_i.T @ S_ov
    V = S_vo @ X_j
    T2 = (U * V).sum()
    return T1 - T2   # CP2K convention: Σc²=1, no factor of 2

if __name__ == "__main__":
    # --- Loop over 19 pairs ---
    print(f"{'Pair':>7}   {'S_11':>10}   {'|S_11|':>8}   {'1 - |S|':>10}")
    print("-" * 50)
    results = []
    for i in range(1, 20):
        j = i + 1
        S_MO = load_swf(f'pairs/frame_{i}_{j}/s_wf_{i}_{j}')
        ci_i = parse_ci_state(f'run_{i}_full.log', state_num=1)
        ci_j = parse_ci_state(f'run_{j}_full.log', state_num=1)
        X_i = build_X(ci_i)
        X_j = build_X(ci_j)
        S = state_overlap(X_i, X_j, S_MO, singlet=True)
        results.append(S)
        print(f"  {i:>2}->{j:<2}   {S:>10.6f}   {abs(S):>8.6f}   {1-abs(S):>10.2e}")
    
    print()
    print("Interpretation:")
    print("  |S_11| ~ 1.000   -> state 1 tracks to state 1 (no crossing)")
    print("  |S_11| < 0.99    -> state 1 has changed character or crossed")
