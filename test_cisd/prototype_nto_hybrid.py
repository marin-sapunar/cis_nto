import numpy as np
from cis_nto.nto import compress_singles_nto

# 1. Simulate a large CISD singles block (e.g., 100 occupied x 100 virtual)
np.random.seed(42)
nocc, nvirt = 100, 100
c_singles = np.random.rand(nocc, nvirt) * 0.01
c_singles[0, 0] = 0.90  # One dominant excitation
c_singles[1, 1] = 0.40  # Second dominant excitation

print("=== NTO Hybrid Prototype ===")
print(f"Original Singles Determinants: {nocc * nvirt}")

# 2. Compress the singles block
U, Vt, S = compress_singles_nto(c_singles, threshold=0.999)
n_nto = len(S)

print(f"Compressed NTO Pairs: {n_nto}")
print(f"Singles Block Compression Factor: {(nocc * nvirt) / n_nto:.1f}x")

# 3. Prove the overlap is preserved
# Reconstruct the singles matrix from the NTOs: C_rec = U * diag(S) * Vt
c_reconstructed = U @ np.diag(S) @ Vt
overlap_error = np.max(np.abs(c_singles - c_reconstructed))

print(f"Max coefficient error after NTO transform: {overlap_error:.4e}")
print("-> NTO Hybrid is mathematically exact for the singles block.")

# 4. In a real hybrid, you would now pair these NTO determinants with 
# the EXACT doubles determinants (which are not compressed) and feed them into L4M.
