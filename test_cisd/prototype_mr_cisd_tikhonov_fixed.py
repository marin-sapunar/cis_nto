import numpy as np
np.set_printoptions(precision=6, suppress=True)

print("=== MR-CISD with Tikhonov-Regularized Schur Update (FIXED) ===\n")

# Multi-reference: 2-det, c0=0.8, c1=0.6
c0, c1 = 0.8, 0.6

S_mo = np.array([
    [0.95, 0.02, 0.01, 0.00],
    [0.02, 0.90, 0.01, 0.01],
    [0.01, 0.01, 0.88, 0.02],
    [0.00, 0.01, 0.02, 0.85],
])

Gamma = np.zeros((2, 2))
Gamma[0, 0] = c0**2 * 2.0
Gamma[1, 1] = c1**2 * 2.0
Gamma[0, 1] = c0 * c1 * 2.0
Gamma[1, 0] = c0 * c1 * 2.0

A_ref = S_mo[:2, :2] @ Gamma @ S_mo[:2, :2].T
print("Gamma:")
print(Gamma)
print(f"det(Gamma) = {np.linalg.det(Gamma):.3e}  (rank-deficient)")
print("\nA_ref:")
print(A_ref)
print(f"det(A_ref) = {np.linalg.det(A_ref):.3e}  (rank-deficient)\n")

# Compute actual rank
n = A_ref.shape[0]
rank = np.linalg.matrix_rank(A_ref, tol=1e-10)
print(f"Matrix size n = {n}, numerical rank r = {rank}")
print(f"Pseudo-determinant divisor: eps^(n - r) = eps^{n - rank}\n")

# Pseudo-determinant: divide by eps^(n - r) = eps^(n - rank)
print("Tikhonov regularization (correct normalization):")
print(f"{'eps':>10} | {'det(A_ref+epsI)':>18} | {'det/eps^(n-r)':>20} | {'converged?':>12}")
print("-" * 70)

pseudo_det_est = None
for eps in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]:
    A_reg = A_ref + eps * np.eye(n)
    det_reg = np.linalg.det(A_reg)
    # Correct divisor: eps^(n - rank)
    ratio = det_reg / (eps ** (n - rank))
    converged = "yes" if pseudo_det_est is not None and abs(ratio - pseudo_det_est) < 0.01 * abs(pseudo_det_est) else "..."
    print(f"{eps:>10.0e} | {det_reg:>18.6e} | {ratio:>20.8f} | {converged:>12}")
    if pseudo_det_est is None or eps < 1e-8:
        pseudo_det_est = ratio

print(f"\nConverged pseudo-determinant = {pseudo_det_est:.8f}")
print("True non-zero eigenvalue of A_ref = ", np.linalg.eigvalsh(A_ref).max())
print("\n-> The ratio converges to the pseudo-determinant, as expected.")
print("-> This is the correct replacement for det(A_ref) in MR-CISD Schur updates.")
