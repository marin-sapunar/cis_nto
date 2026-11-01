import numpy as np
import time
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score

def load_amyloid_case(filename):
    with open(filename, 'r') as f:
        header = f.readline().split()
        n_mo, n_occ, n_dets = int(header[0]), int(header[1]), int(header[2])
        S_mo = np.zeros((n_mo, n_mo))
        for i in range(n_mo):
            S_mo[i, :] = [float(x) for x in f.readline().split()]
        dets_a = np.zeros((n_occ, n_dets), dtype=int)
        dets_b = np.zeros((n_occ, n_dets), dtype=int)
        coeffs = np.zeros(n_dets)
        for k in range(n_dets):
            dets_a[:, k] = [int(x)-1 for x in f.readline().split()]
            dets_b[:, k] = [int(x)-1 for x in f.readline().split()]
            coeffs[k] = float(f.readline().split()[0])
    return S_mo, n_occ, n_dets, dets_a, dets_b, coeffs

def compute_det_overlaps(S_mo, n_occ, dets_a):
    # For closed-shell amyloid data, alpha and beta blocks are identical.
    # det_{ij} = det(S_mo[occ_i, occ_j])^2
    n_dets = dets_a.shape[1]
    D = np.zeros((n_dets, n_dets))
    for I in range(n_dets):
        for J in range(n_dets):
            M = S_mo[np.ix_(dets_a[:, I], dets_a[:, J])]
            D[I, J] = np.linalg.det(M) ** 2
    return D

def contributions_from_D(D, c):
    # S_full = c^T D c
    # S_drop_i = S_full - 2 c_i (Dc)_i + c_i^2 D_ii
    S_full = c @ D @ c
    Dc = D @ c
    contrib = np.zeros(len(c))
    for i in range(len(c)):
        S_drop = S_full - 2*c[i]*Dc[i] + c[i]**2 * D[i,i]
        contrib[i] = abs(S_full - S_drop)
    return S_full, contrib

# ---- Main ----
filename = 'amyloid_t3_d0.05_n50.txt'
S_mo, n_occ, n_dets, dets_a, dets_b, coeffs = load_amyloid_case(filename)
print(f"Loaded {filename}: {n_dets} determinants, {n_occ} occ, {S_mo.shape[0]} MOs")

t0 = time.time()
D = compute_det_overlaps(S_mo, n_occ, dets_a)
print(f"Computed {n_dets}x{n_dets} det overlaps in {time.time()-t0:.1f}s")

S_full, contrib = contributions_from_D(D, coeffs)
print(f"S_full = {S_full:.10f}\n")

# ---- Features ----
features = []
for i in range(n_dets):
    c_i = abs(coeffs[i])
    n_exc = (dets_a[:, i] != np.arange(n_occ)).sum()
    occ = dets_a[:, i]
    diag_ovlp = S_mo[np.ix_(occ, occ)]
    diag_sum = np.sum(diag_ovlp**2)
    max_ovlp = np.max(np.abs(diag_ovlp))
    features.append([c_i, c_i**2, n_exc, diag_sum, max_ovlp])
features = np.array(features)

from numpy import corrcoef
corr = corrcoef(features[:, 1], contrib)[0,1]
print(f"Correlation |c|^2 vs true contribution: {corr:.3f}")
print("(If < 1, coefficient-based truncation is suboptimal)\n")

# ---- Train ----
model = GradientBoostingRegressor(n_estimators=200, max_depth=3, random_state=42)
scores = cross_val_score(model, features, contrib, cv=5, scoring='r2')
print(f"5-fold CV R^2: {scores.mean():.3f} ± {scores.std():.3f}")
model.fit(features, contrib)

# ---- Compare truncation strategies ----
print(f"\n{'K':>4} | {'Baseline |c|^2 err':>20} | {'ML-guided err':>20} | {'Improvement':>12}")
print("-" * 72)
for K in [5, 10, 15, 20, 30, 40]:
    idx_b = np.argsort(features[:, 1])[::-1][:K]
    S_b = coeffs[idx_b] @ D[np.ix_(idx_b, idx_b)] @ coeffs[idx_b]
    err_b = abs(S_full - S_b)

    pred = model.predict(features)
    idx_ml = np.argsort(pred)[::-1][:K]
    S_ml = coeffs[idx_ml] @ D[np.ix_(idx_ml, idx_ml)] @ coeffs[idx_ml]
    err_ml = abs(S_full - S_ml)

    imp = (err_b - err_ml) / err_b * 100 if err_b > 1e-15 else 0
    print(f"{K:>4} | {err_b:>20.4e} | {err_ml:>20.4e} | {imp:>10.1f}%")

print("\n-> If ML-guided error < baseline error at same K,")
print("   ML successfully identifies which determinants matter for the overlap.")
