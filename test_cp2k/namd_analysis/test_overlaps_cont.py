import numpy as np, pickle
from state_overlap import build_X, state_overlap
from read_molden_mo import read_molden_mo

def lowdin_inv_sqrt(S):
    vals, vecs = np.linalg.eigh(S)
    return vecs @ np.diag(1.0 / np.sqrt(np.maximum(vals, 1e-12))) @ vecs.T

for i in 1, 2, 3:
    j = i + 1
    S_AO = np.loadtxt(f's_ao_cont_{i}_{j}', skiprows=1)
    C_A = read_molden_mo(f'p4-di-PBE_cont_{i}-MOS-1_0.molden')
    C_B = read_molden_mo(f'p4-di-PBE_cont_{j}-MOS-1_0.molden')
    
    X_AB = C_A.T @ S_AO @ C_B
    X_lowdin = lowdin_inv_sqrt(C_A.T @ S_AO @ C_A) @ X_AB @ lowdin_inv_sqrt(C_B.T @ S_AO @ C_B)
    
    with open(f'states_cont_{i}.pkl', 'rb') as f: states_A = pickle.load(f)
    with open(f'states_cont_{j}.pkl', 'rb') as f: states_B = pickle.load(f)
        
    S_11 = state_overlap(build_X(states_A[1]['ci']), build_X(states_B[1]['ci']), X_lowdin, singlet=False)
    sig = np.linalg.svd(X_lowdin, compute_uv=False)
    
    print(f"--- Pair {i} -> {j} ---")
    print(f"Sigma Max: {sig.max():.8f}  |  S_11: {abs(S_11):.4f}")
