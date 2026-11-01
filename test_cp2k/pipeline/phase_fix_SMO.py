import numpy as np
import pickle

def phase_fix_S_MO(S_MO):
    """Flip MO signs so S_MO[i,i] > 0 for all i (CP2K molden phase ambiguity)."""
    S = S_MO.copy()
    n_flips = 0
    for i in range(S.shape[0]):
        if S[i, i] < 0:
            S[i, :] *= -1
            S[:, i] *= -1
            n_flips += 1
    return S, n_flips

# Phase-fix all 19 S_MO files in place (keep originals)
import os
for i in range(1, 20):
    j = i + 1
    f = f'pairs/frame_{i}_{j}/s_wf_{i}_{j}'
    S = np.loadtxt(f, skiprows=1)
    S_fixed, n = phase_fix_S_MO(S)
    with open(f'{f}_phasefixed', 'w') as fh:
        fh.write(f'{S_fixed.shape[0]}\n')
        np.savetxt(fh, S_fixed, fmt='%.10e')
    print(f'{i:>2}→{j:<2}: flipped {n:>3d} of 205 MOs')
