import numpy as np

# 1. Load the overlap matrix
S = np.array([
    [ 0.987557,  0.000000,  0.000000],
    [-0.000000,  0.545017, -0.002047],
    [ 0.000000, -0.002478, -0.711793]
])

dt_au = 0.5 * 41.341  # 0.5 fs in atomic units

# 2. Compute TDCs
n_states = S.shape[0]
TDC = np.zeros((n_states, n_states))
for i in range(n_states):
    for j in range(n_states):
        if i != j:
            TDC[i, j] = (S[i, j] - S[j, i]) / (2.0 * dt_au)

# 3. Simulate a realistic superposition state
# We are in a mixed state of S1 and S2 (coherence is built up)
c = np.zeros(n_states)
c[1] = np.sqrt(0.5)
c[2] = np.sqrt(0.5)

# 4. Compute Hopping Probability S1 -> S2
current_state = 1
target_state = 2

prob = 2.0 * np.real(c[current_state] * np.conj(c[target_state]) * TDC[current_state, target_state]) * dt_au

print(f"FSSH Hopping Probability S1 -> S2: {prob:.6e}")

if prob > 0:
    print("-> Hop is ALLOWED (positive probability)")
    print(f"-> The overlap matrix successfully generates a hopping probability of {prob*100:.4f}% per step.")
else:
    print("-> Hop is forbidden")
