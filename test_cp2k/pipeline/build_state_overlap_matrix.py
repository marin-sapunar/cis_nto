import pickle
import numpy as np
from state_overlap import load_swf, build_X, state_overlap

with open('all_states_pairphase.pkl', 'rb') as f:
    all_states = pickle.load(f)

for i in range(1, 20):
    j = i + 1
    S_MO = load_swf(f'pairs/frame_{i}_{j}/s_wf_{i}_{j}_phasefixed')
    states_i = all_states[i]
    states_j = all_states[j]

    n_states = min(len(states_i), len(states_j), 10)
    S_IJ = np.zeros((n_states, n_states))

    for I in range(1, n_states + 1):
        X_I = build_X(states_i[I]['ci'])
        for J in range(1, n_states + 1):
            X_J = build_X(states_j[J]['ci'])
            S_IJ[I-1, J-1] = state_overlap(X_I, X_J, S_MO, singlet=False)

    np.save(f'pairs/frame_{i}_{j}/S_IJ_{i}_{j}.npy', S_IJ)

    # Check that diagonal is ~1
    diag = np.diag(S_IJ)
    off_diag_max = np.abs(S_IJ - np.diag(diag)).max()
    print(f"{i:>2}->{j:<2}:  diag ∈ [{diag.min():.4f}, {diag.max():.4f}]  "
          f"max off-diag = {off_diag_max:.4f}")
