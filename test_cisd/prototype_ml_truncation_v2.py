import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

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

import glob
files = sorted(glob.glob('amyloid_t*_n*.txt'))
all_X, all_y = [], []
for f in files:
    S_mo, n_occ, dets_a, coeffs = load_case(f)
    n_dets = len(coeffs)
    # Det overlaps in blocks
    D = np.zeros((n_dets, n_dets))
    for I in range(n_dets):
        for J in range(n_dets):
            D[I, J] = np.linalg.det(S_mo[np.ix_(dets_a[:, I], dets_a[:, J])])**2
    S_full = coeffs @ D @ coeffs
    Dc = D @ coeffs
    for i in range(n_dets):
        S_drop = S_full - 2*coeffs[i]*Dc[i] + coeffs[i]**2 * D[i,i]
        contrib = abs(S_full - S_drop)
        occ = dets_a[:, i]
        diag = S_mo[np.ix_(occ, occ)]
        n_exc = (occ != np.arange(n_occ)).sum()
        all_X.append([abs(coeffs[i]), coeffs[i]**2, n_exc, np.sum(diag**2), np.max(np.abs(diag))])
        all_y.append(contrib)

X = np.array(all_X)
y = np.array(all_y)
print(f"Total samples across all cases: {len(X)}")
print(f"Correlation |c|^2 vs contribution: {np.corrcoef(X[:,1], y)[0,1]:.3f}")

pipe = Pipeline([('scale', StandardScaler()), ('ridge', RidgeCV(alphas=[0.1, 1, 10]))])
scores = cross_val_score(pipe, X, y, cv=5, scoring='r2')
print(f"5-fold CV R^2 (Ridge): {scores.mean():.3f} +/- {scores.std():.3f}")
