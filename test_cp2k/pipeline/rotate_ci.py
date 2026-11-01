"""
Rotate a CI vector between two frames using the pseudo-inverse
biorthogonal transform. Norm-preserving by construction.

Usage:
    python rotate_ci.py <S_MO_file> <ci_log> [n_occ] [n_virt]
"""
import sys
import re
import numpy as np

def load_swf(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        return np.loadtxt(f, max_rows=n)

def parse_ci_state1(log_file):
    """Return list of (occ, virt, amplitude) for state 1."""
    with open(log_file) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if "Excitation analysis" in l)
    state_re = re.compile(r'^\s+1\s+\d+\.\d+ eV')
    next_re = re.compile(r'^\s+\d+\s+\d+\.\d+ eV')
    i = start
    while not state_re.match(lines[i]):
        i += 1
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

def build_X(ci, n_occ, n_virt):
    X = np.zeros((n_occ, n_virt))
    for occ, virt, amp in ci:
        v = virt - (n_occ + 1)
        if 0 <= v < n_virt:
            X[occ-1, v] = amp
    return X

def rotate_ci(X, S_MO, n_occ):
    """Biorthogonal rotation using pseudo-inverse."""
    S_occ = S_MO[:n_occ, :n_occ]
    S_virt = S_MO[n_occ:, n_occ:]
    return np.linalg.pinv(S_occ) @ X @ np.linalg.pinv(S_virt).T

if __name__ == "__main__":
    s_mo_file = sys.argv[1]
    log_file = sys.argv[2]
    n_occ = int(sys.argv[3]) if len(sys.argv) > 3 else 155
    n_virt = int(sys.argv[4]) if len(sys.argv) > 4 else 50

    S_MO = load_swf(s_mo_file)
    ci = parse_ci_state1(log_file)
    X = build_X(ci, n_occ, n_virt)
    X_t = rotate_ci(X, S_MO, n_occ)

    norm_o = np.linalg.norm(X)
    norm_t = np.linalg.norm(X_t)
    print(f"Frame pair: {s_mo_file}")
    print(f"  ||X_orig||  = {norm_o:.8f}")
    print(f"  ||X_trans|| = {norm_t:.8f}")
    print(f"  ratio       = {norm_t/norm_o:.8f}")
    print(f"  error       = {abs(norm_t/norm_o - 1):.2e}")
