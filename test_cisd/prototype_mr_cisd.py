import numpy as np

# 1. Define a simple 4x4 MO overlap matrix between t1 and t2
# (This represents the MO overlaps of the active space)
S_mo = np.array([
    [ 0.95,  0.02,  0.01,  0.00],
    [ 0.02,  0.90,  0.01,  0.01],
    [ 0.01,  0.01,  0.88,  0.02],
    [ 0.00,  0.01,  0.02,  0.85]
])

# 2. Single-reference (Standard CISD): 
# The reference is just the HF determinant (first 2 orbitals occupied).
# The reference matrix is the top-left 2x2 block of S_mo.
A_ref_sr = S_mo[:2, :2]
det_A_sr = np.linalg.det(A_ref_sr)

print("=== MR-CISD Reference Matrix Prototype ===")
print("Single-Reference (HF) A_ref determinant:", det_A_sr)

# 3. Multi-reference (MR-CISD):
# The reference is a linear combination of the HF determinant and a doubly-excited determinant.
# The multi-reference density matrix is: Gamma = sum_IJ c_I c_J * <Phi_I|E|Phi_J>
# For a 2-determinant mix, Gamma = c0^2 * Gamma_0 + c1^2 * Gamma_1 + 2*c0*c1*Gamma_01

# Let's say the active space has 2 electrons in 2 orbitals.
# |Phi_0> = |1^2>
# |Phi_1> = |2^2>
# The density matrix in the active space is:
c0, c1 = 0.8, 0.6  # Coefficients from CASSCF

# Active space density matrix (in the 2 active orbitals)
gamma_act = np.zeros((2, 2))
# |Phi_0> contributes to orbital 1
gamma_act[0, 0] += c0**2 * 2.0
# |Phi_1> contributes to orbital 2
gamma_act[1, 1] += c1**2 * 2.0
# Off-diagonal coherence terms
gamma_act[0, 1] += 2 * c0 * c1 * 1.0
gamma_act[1, 0] += 2 * c0 * c1 * 1.0

print("\nActive Space Density Matrix (Gamma):")
print(gamma_act)

# The multi-reference A_ref matrix is built by contracting S_mo with Gamma
# A_ref_ij = sum_kl S_mo_ik * Gamma_kl * S_mo_jl
A_ref_mr = S_mo[:2, :2] @ gamma_act @ S_mo[:2, :2].T

print("\nMulti-Reference A_ref matrix:")
print(A_ref_mr)
print("Multi-Reference determinant:", np.linalg.det(A_ref_mr))

print("\n-> The reference matrix is no longer a simple MO overlap block; it includes density correlations from the active space.")
