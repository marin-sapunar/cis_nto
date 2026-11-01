import numpy as np
import pickle
from state_overlap import build_X, state_overlap
from read_molden_mo import read_molden_mo

def lowdin_inv_sqrt(S):
    vals, vecs = np.linalg.eigh(S)
    return vecs @ np.diag(1.0 / np.sqrt(np.maximum(vals, 1e-12))) @ vecs.T

# NAMD Parameters
dt = 1.0 # fs (since frames are consecutive 1 fs apart)
print(f"{'Time (fs)':<10} | {'State':<6} | {'Energy (eV)':<12} | {'Max NAC (fs⁻¹)':<15} | {'Max Hop Prob':<12}")
print("-" * 75)

for i in range(1, 20):
    j = i + 1
    
    # 1. Load data for pair i -> j
    S_AO = np.loadtxt(f's_ao_cont_{i}_{j}', skiprows=1)
    C_A = read_molden_mo(f'p4-di-PBE_cont_{i}-MOS-1_0.molden')
    C_B = read_molden_mo(f'p4-di-PBE_cont_{j}-MOS-1_0.molden')
    
    # 2. Build cross-overlap S_MO
    X_AB = C_A.T @ S_AO @ C_B
    X_lowdin = lowdin_inv_sqrt(C_A.T @ S_AO @ C_A) @ X_AB @ lowdin_inv_sqrt(C_B.T @ S_AO @ C_B)
    
    # 3. Load CI vectors
    with open(f'states_cont_{i}.pkl', 'rb') as f: states_A = pickle.load(f)
    with open(f'states_cont_{j}.pkl', 'rb') as f: states_B = pickle.load(f)
    
    # 4. Build the 15x15 S_IJ matrix for this time step
    n_states = 15
    S_IJ = np.zeros((n_states, n_states))
    for I in range(1, n_states + 1):
        X_I = build_X(states_A[I]['ci'])
        for J in range(1, n_states + 1):
            X_J = build_X(states_B[J]['ci'])
            S_IJ[I-1, J-1] = state_overlap(X_I, X_J, X_lowdin, singlet=False)
            
    # 5. Compute Time-Derivative NAC (off-diagonals)
    # d_IJ = (S_IJ - S_JI) / (2 * dt)
    NAC = (S_IJ - S_IJ.T) / (2.0 * dt)
    np.fill_diagonal(NAC, 0.0)
    
    # 6. Compute Hopping Probabilities (assuming a generic velocity v=1 for this test)
    # P = max(0, -2 * Re(d_IJ * v) * dt)
    # For this test, we just want to see if the probabilities are physically reasonable.
    Hop_Prob = np.maximum(0.0, -2.0 * NAC * 1.0 * dt)
    
    # Print max NAC and Hop Prob for the dominant state (State 1)
    # Since we track State 1 for this test:
    max_nac = np.max(np.abs(NAC[0, 1:]))
    max_hop = np.max(Hop_Prob[0, 1:])
    energy = states_A[1]['energy']
    
    print(f"{(i-1)*dt:<10.1f} | {'S1':<6} | {energy:<12.4f} | {max_nac:<15.4f} | {max_hop:<12.4f}")

