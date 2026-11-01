import pickle, numpy as np
from scipy.optimize import linear_sum_assignment
from state_overlap import load_swf, build_X, state_overlap

def orthog_lowdin(S):
    U, sv, Vt = np.linalg.svd(S, full_matrices=False)
    return U @ Vt, sv

def phase_match_rotation(S):
    S = S.copy(); n = S.shape[0]
    r, c = linear_sum_assignment(-(S**2)); cm = np.empty(n, dtype=int); cm[r] = c
    sw = np.zeros(n, dtype=bool)
    for st in range(n):
        if S[st, cm[st]] < 0:
            S[:, cm[st]] *= -1; sw[cm[st]] = ~sw[cm[st]]
    vis = np.zeros(n, dtype=bool)
    for start in range(n):
        if vis[start]: continue
        lp = []; cur = start
        while not vis[cur]:
            vis[cur] = True; lp.append(cur); cur = cm[cur]
        if len(lp) > 1 and np.linalg.det(S[np.ix_(lp, lp)]) < 0:
            S[:, lp[0]] *= -1; sw[lp[0]] = ~sw[lp[0]]
    return S, cm, sw

with open('all_states_pairphase.pkl','rb') as f: all_states = pickle.load(f)
print(f"{'pair':<8} {'raw':>10} {'oo':>10} {'oo+mp':>10} {'swaps':>5}")
for i in range(1, 20):
    j = i+1
    try:
        S_MO = load_swf(f'pairs/frame_{i}_{j}/s_wf_{i}_{j}_phasefixed')
    except FileNotFoundError:
        print(f"{i}->{j:<6} MISSING"); continue
    si, sj = all_states[i], all_states[j]
    n = min(len(si), len(sj), 10)
    S = np.zeros((n, n))
    for I in range(1, n+1):
        XI = build_X(si[I]['ci'])
        for J in range(1, n+1):
            XJ = build_X(sj[J]['ci'])
            S[I-1,J-1] = state_overlap(XI, XJ, S_MO, singlet=False)
    raw_off = np.abs(S - np.diag(np.diag(S))).max()
    Q, _ = orthog_lowdin(S)
    oo_off  = np.abs(Q - np.diag(np.diag(Q))).max()
    Mp, _, sw = phase_match_rotation(Q)
    mp_off  = np.abs(Mp - np.diag(np.diag(Mp))).max()
    print(f"{i}->{j:<6} {raw_off:>10.4e} {oo_off:>10.4e} {mp_off:>10.4e} {int(sw.sum()):>5}")
