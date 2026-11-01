import numpy as np

n_states = 3
n_steps = 50
dt = 0.5 # fs

print("Generating Ethylene Conical Intersection Trajectory...")

for step in range(n_steps):
    # Simulate the C=C twist angle (0 to 90 degrees)
    theta = (step / n_steps) * (np.pi / 2)
    
    # At the CI, State 2 (bright) mixes with State 1 (ground) and State 3 (dark)
    # The overlap matrix elements change as a function of the twist angle
    S = np.eye(n_states)
    
    # Off-diagonal couplings grow as we approach the CI (theta = 45 deg)
    coupling_12 = 0.05 * np.sin(2 * theta) 
    coupling_13 = 0.01 * np.sin(2 * theta)
    coupling_23 = 0.02 * np.sin(2 * theta)
    
    # Construct a NON-SYMMETRIC overlap matrix (S_ij != S_ji) to simulate time evolution
    S[0, 1] = coupling_12 + 0.01
    S[1, 0] = coupling_12 - 0.01
    S[0, 2] = coupling_13 + 0.005
    S[2, 0] = coupling_13 - 0.005
    S[1, 2] = coupling_23 + 0.002
    S[2, 1] = coupling_23 - 0.002
    
    # Add a tiny decay to the diagonal (nuclear motion)
    for i in range(n_states):
        S[i, i] = 1.0 - (step * 0.0005)
        
    np.savetxt(f"st_adi_{step}.txt", S, fmt='%.8f')

print("✅ Generated 50 overlap matrices for Ethylene CI.")
