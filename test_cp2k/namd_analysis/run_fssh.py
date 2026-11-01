import numpy as np
import pickle
from state_overlap import build_X, state_overlap
from read_molden_mo import read_molden_mo

def lowdin_inv_sqrt(S):
    vals, vecs = np.linalg.eigh(S)
    return vecs @ np.diag(1.0 / np.sqrt(np.maximum(vals, 1e-12))) @ vecs.T

# FSSH simulation
current_state = 1
dt = 1.0
np.random.seed(42) # For reproducibility

print(f"Starting FSSH in State {current_state}...")
print(f"{'Time (fs)':<10} | {'Current State':<13} | {'Energy (eV)':<12} | {'Hopped To?'}")
print("-" * 60)

for i in range(1, 20):
    j = i + 1
    S_AO = np.loadtxt(f's_ao_cont_{i}_{j}', skiprows=1)
    C_A = read_molden_mo(f'p4-di-PBE_cont_{i}-MOS-1_0.molden')
    C_B = read_molden_mo(f'p4-di-PBE_cont_{j}-MOS-1_0.molden')
    
    X_AB = C_A.T @ S_AO @ C_B
    X_lowdin = lowdin_inv_sqrt(C_A.T @ S_AO @ C_A) @ X_AB @ lowdin_inv_sqrt(C_B.T @ S_AO @ C_B)
    
    with open(f'states_cont_{i}.pkl', 'rb') as f: states_A = pickle.load(f)
    with open(f'states_cont_{j}.pkl', 'rb') as f: states_B = pickle.load(f)
    
    # Calculate probabilities for the CURRENT state we are in
    X_I = build_X(states_A[current_state]['ci'])
    probs = []
    for J in range(1, 16):
        X_J = build_X(states_B[J]['ci'])
        S = state_overlap(X_I, X_J, X_lowdin, singlet=False)
        if J != current_state:
            # Simplified FSSH probability for this validation test
            prob = max(0.0, -2.0 * S * 1.0 * dt) 
            probs.append((J, prob))
            
    # Check if we hop
    hopped = False
    for J, prob in probs:
        if np.random.rand() < prob:
            print(f"{(i-1)*dt:<10.1f} | {current_state:<13} | {states_A[current_state]['energy']:<12.4f} | HOPPED TO S{J}!")
            current_state = J
            hopped = True
            break
            
    if not hopped:
         print(f"{(i-1)*dt:<10.1f} | {current_state:<13} | {states_A[current_state]['energy']:<12.4f} | No hop")

