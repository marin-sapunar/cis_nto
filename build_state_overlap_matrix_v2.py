import pickle
import numpy as np
from scipy.optimize import linear_sum_assignment
from state_overlap import load_swf, build_X, state_overlap

ORTH_OVERLAP = False
MATCH_PHASE  = False
RANK_TOL     = 1e-10

def orthog_lowdin(S):
    U, sv, Vt = np.linalg.svd(S, full_matrices=False)
    return U @ Vt, sv

def phase_match_assigned_rotation(S):
    S = S.copy()
    n = S.shape[0]
    row_ind, col_ind = linear_sum_assignment(-(S ** 2))
    cmatch = np.empty(n, dtype=int)
    cmatch[row_ind] = col_ind
    swap = np.zeros(n, dtype=bool)
    for st in range(n):
        if S[st, cmatch[st]] < 0:
            S[:, cmatch[st]] *= -1
            swap[cmatch[st]] = ~swap[cmatch[st]]
    visited = np.zeros(n, dtype=bool)
    for start in range(n):
        if visited[start]:
            continue
        loop = []
        cur = start
        while not visited[cur]:
            visited[cur] = True
            loop.append(cur)
            cur = cmatch[cur]
        if len(loop) <= 1:
            continue
        if np.linalg.det(S[np.ix_(loop, loop)]) < 0:
            S[:, loop[0]] *= -1
            swap[loop[0]] = ~swap[loop[0]]
    return S, cmatch, swap

with open('all_states_pairphase.pkl', 'rb') as f:
    all_states = pickle.load(f)

for i in range(1, 20):
    j = i + 1
    S_MO = load_swf(f'pairs/frame_{i}_{j}/s_wf_{i}_{j}_phasefixed')
    states_i, states_j = all_states[i], all_states[j]
    n_states = min(len(states_i), len(states_j), 10)
    S_IJ = np.zeros((n_states, n_states))
    for I in range(1, n_states + 1):
        X_I = build_X(states_i[I]['ci'])
        for J in range(1, n_states + 1):
            X_J = build_X(states_j[J]['ci'])
            S_IJ[I-1, J-1] = state_overlap(X_I, X_J, S_MO, singlet=False)
    sv = None
    rank = n_states
    if ORTH_OVERLAP:
        S_IJ, sv = orthog_lowdin(S_IJ)
        rank = int(np.sum(sv > RANK_TOL * sv[0])) if sv.size and sv[0] > 0 else 0
    n_swaps = 0
    if MATCH_PHASE:
        S_IJ, cmatch, swap = phase_match_assigned_rotation(S_IJ)
        n_swaps = int(swap.sum())
    np.save(f'pairs/frame_{i}_{j}/S_IJ_{i}_{j}.npy', S_IJ)
    diag = np.abs(np.diag(S_IJ))
    off  = np.abs(S_IJ - np.diag(np.diag(S_IJ))).max()
    extra = ""
    if ORTH_OVERLAP:
        orth_err = np.linalg.norm(S_IJ.T @ S_IJ - np.eye(n_states), ord='fro')
        extra += f"  orth_err={orth_err:.1e} rank={rank}"
    if MATCH_PHASE:
        extra += f"  n_swaps={n_swaps}"
    print(f"{i:>2}->{j:<2}: diag ∈ [{diag.min():.6f},{diag.max():.6f}]  "
          f"max off = {off:.4e}{extra}")

# Trigger PR diff refresh

# Trigger PR diff refresh
