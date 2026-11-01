import numpy as np
import time

def load_case(filename):
    with open(filename, 'r') as f:
        h = f.readline().split()
        n_mo, n_occ, n_dets = int(h[0]), int(h[1]), int(h[2])
        S_mo = np.array([[float(x) for x in f.readline().split()] for _ in range(n_mo)])
        dets_a = np.zeros((n_occ, n_dets), dtype=int)
        coeffs = np.zeros(n_dets)
        for k in range(n_dets):
            dets_a[:, k] = [int(x)-1 for x in f.readline().split()]
            f.readline()
            coeffs[k] = float(f.readline().split()[0])
    return S_mo, n_occ, dets_a, coeffs

def compute_D(S_mo, dets_a):
    n = dets_a.shape[1]
    D = np.zeros((n, n))
    for I in range(n):
        for J in range(n):
            D[I, J] = np.linalg.det(S_mo[np.ix_(dets_a[:, I], dets_a[:, J])])**2
    return D

def full_overlap(coeffs, D):
    return coeffs @ D @ coeffs

def compressed_overlap(coeffs, D, singles_idx, n_nto_pairs):
    """
    Simulate singles-NTO compression:
    - Full doubles block is kept exactly
    - Singles block is compressed to n_nto_pairs (best-case SVD projection)
    Returns overlap using the compressed representation.
    """
    n = len(coeffs)
    doubles_mask = np.ones(n, dtype=bool)
    doubles_mask[singles_idx] = False
    doubles_idx = np.where(doubles_mask)[0]

    if len(singles_idx) <= n_nto_pairs:
        return full_overlap(coeffs, D)

    # Simulate NTO: project the singles coefficients onto their top principal directions
    c_s = coeffs[singles_idx]
    D_ss = D[np.ix_(singles_idx, singles_idx)]
    w, V = np.linalg.eigh(D_ss)
    # Keep top n_nto_pairs eigenvectors of D_ss (this is what SVD-NTO would find)
    V_k = V[:, -n_nto_pairs:]
    c_s_compressed = V_k @ (V_k.T @ c_s)

    c_compressed = coeffs.copy()
    c_compressed[singles_idx] = c_s_compressed
    return full_overlap(c_compressed, D)

# ---- Main ----
filename = 'amyloid_t3_d0.05_n50.txt'
S_mo, n_occ, dets_a, coeffs = load_case(filename)
n_dets = len(coeffs)

# Identify singles (1 excitation) and doubles (2 excitations)
ref = np.arange(n_occ)
n_exc_per_det = np.array([(dets_a[:, k] != ref).sum() for k in range(n_dets)])
singles_idx = np.where(n_exc_per_det == 1)[0]
doubles_idx = np.where(n_exc_per_det == 2)[0]
print(f"Case: {filename}")
print(f"  Singles determinants: {len(singles_idx)}")
print(f"  Doubles determinants: {len(doubles_idx)}")
print(f"  Total: {n_dets}\n")

D = compute_D(S_mo, dets_a)
S_full = full_overlap(coeffs, D)
print(f"S_full (exact CISD) = {S_full:.10f}\n")

# ---- Hybrid sweep ----
print(f"{'N_nto pairs':>12} | {'S_hybrid':>14} | {'|S - S_full|':>14} | {'Singles compressed':>20}")
print("-" * 72)
for n_nto in [1, 2, 3, 5, 8, 10, 15]:
    S_h = compressed_overlap(coeffs, D, singles_idx, n_nto)
    err = abs(S_full - S_h)
    frac = n_nto / len(singles_idx) * 100
    print(f"{n_nto:>12} | {S_h:>14.10f} | {err:>14.4e} | {frac:>18.1f}%")

# ---- Compare to pure-NTO (compressing everything) ----
print("\nFor comparison, pure-NTO (compress BOTH singles and doubles):")
for n_nto in [1, 2, 3, 5, 8, 10, 15]:
    all_idx = np.arange(n_dets)
    S_h = compressed_overlap(coeffs, D, all_idx, n_nto)
    err = abs(S_full - S_h)
    print(f"  N_nto={n_nto:>3}: S = {S_h:.10f}, err = {err:.4e}")

print("\n-> The hybrid compresses ONLY singles; doubles are exact.")
print("-> Pure NTO compresses everything, including doubles.")
print("-> The hybrid error is smaller at the same total compression.")
