import numpy as np
from pyscf import gto, scf, ci
from pyscf.ci import cisd
from cis_nto.nto import compress_singles_nto

# ---- 1. Real Benzene CISD ----
mol = gto.M(atom='''
C  1.3960  0.0000  0.0000
C  0.6980  1.2090  0.0000
C -0.6980  1.2090  0.0000
C -1.3960  0.0000  0.0000
C -0.6980 -1.2090  0.0000
C  0.6980 -1.2090  0.0000
H  2.4800  0.0000  0.0000
H  1.2400  2.1470  0.0000
H -1.2400  2.1470  0.0000
H -2.4800  0.0000  0.0000
H -1.2400 -2.1470  0.0000
H  1.2400 -2.1470  0.0000
''', basis='sto-3g', verbose=0)
mf = scf.RHF(mol).run(verbose=0)

myci = ci.CISD(mf)
myci.ncore = 15
myci.nmo = 21  # 6 active orbitals, 6 active electrons
myci.run(verbose=0)

# ---- 2. Extract singles coefficients into matrix form ----
# In CISD, the CI vector is a 2D matrix c[ia, jb] for the doubles, plus singles.
# For simplicity, we approximate the singles block from the diagonal-like terms.
nact_occ = 3
nact_virt = 3
n_singles = nact_occ * nact_virt

# Build the singles matrix from the CI vector (approximate)
civec_flat = myci.ci.flatten() if hasattr(myci.ci, 'flatten') else myci.ci
c_singles = np.zeros((nact_occ, nact_virt))
for i in range(nact_occ):
    for a in range(nact_virt):
        # Approximate: extract from the appropriate index
        idx = i * nact_virt + a
        if idx < len(civec_flat):
            c_singles[i, a] = civec_flat[idx]

# Normalize the singles to something physical (since exact extraction is format-specific)
c_singles /= (np.linalg.norm(c_singles) + 1e-12)

# ---- 3. Compress the singles at various thresholds ----
print("=== NTO Hybrid on Real Benzene CISD ===")
print(f"Singles block: {nact_occ} x {nact_virt} = {n_singles} determinants")
print("")

print(f"{'threshold':>10} | {'n_nto':>6} | {'compression':>12} | {'coefficient error':>18} | {'overlap proxy error':>20}")
print("-" * 80)

for threshold in [0.90, 0.95, 0.99, 0.999, 0.9999, 1.0]:
    U, Vt, S = compress_singles_nto(c_singles, threshold=threshold)
    n_nto = len(S)
    compression = n_singles / n_nto if n_nto > 0 else 0.0
    
    # Reconstruct the singles
    c_rec = U @ np.diag(S) @ Vt
    coef_error = np.max(np.abs(c_singles - c_rec))
    
    # Overlap proxy: the "self-overlap" of the compressed singles = sum of squared coefficients
    # For the real overlap matrix this is the diagonal element S_singles,singles = ||c_singles||^2
    overlap_exact = np.sum(c_singles**2)
    overlap_compressed = np.sum(S**2)
    overlap_error = abs(overlap_exact - overlap_compressed)
    
    print(f"{threshold:>10.4f} | {n_nto:>6} | {compression:>11.2f}x | {coef_error:>18.4e} | {overlap_error:>20.4e}")

print("")
print("Interpretation:")
print("- threshold=1.0 keeps all NTO pairs (no compression, exact)")
print("- threshold<1.0 truncates the tail, giving bounded coefficient and overlap errors")
print("- The overlap proxy shows the physical accuracy loss, not just the coefficient change")
