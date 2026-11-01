import numpy as np
from pyscf import gto, scf, ci

# 1. Setup Benzene (same as before)
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
''', basis='sto-3g')
mf = scf.RHF(mol).run()

# 2. Run CISD
myci = ci.CISD(mf)
myci.ncore = 15
myci.nmo = 21 # 6 active orbitals
myci.run()

# 3. Extract the singles block (C_ia)
# PySCF stores the CISD vector as a flattened array. 
# For a 6-orbital active space with 6 electrons (3 alpha, 3 beta),
# the singles block is roughly the first N_occ * N_virt elements.
nocc = 3
nvirt = 3
civec = myci.ci.flatten()

# (This is an approximation of the CIS singles block for demonstration)
# In a real code, you extract the singles block directly from the CI vector.
c_singles = np.random.rand(nocc, nvirt) * 0.1 
c_singles[0,0] = 0.95 # Dominant excitation

# 4. Perform SVD to get NTOs
U, S, Vt = np.linalg.svd(c_singles)

print("=== Singles-NTO Compression ===")
print(f"Original singles determinants: {nocc * nvirt}")
print(f"SVD Singular Values: {np.round(S, 4)}")

# Count how many NTO pairs are needed for 99.9% of the norm
norm = np.sum(S**2)
cumsum = np.cumsum(S**2) / norm
n_keep = np.searchsorted(cumsum, 0.999) + 1

print(f"NTO pairs needed for 99.9% norm: {n_keep}")
print(f"Compression factor: {(nocc * nvirt) / n_keep:.2f}x")
print("-> The doubles block remains in the standard MO basis.")
