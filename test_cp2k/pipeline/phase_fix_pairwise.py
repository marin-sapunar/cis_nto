import pickle
import numpy as np

def ci_dot(ci_a, ci_b):
    """Dot product of two sparse CI vectors."""
    da = {(o, v): a for o, v, a in ci_a}
    db = {(o, v): a for o, v, a in ci_b}
    return sum(da[k] * db[k] for k in set(da) & set(db))

with open('all_states.pkl', 'rb') as f:
    all_states = pickle.load(f)

# Pairwise phase fix: match frame N+1's state S to frame N's state S
# Propagate sign forward along the trajectory
for s in range(1, 11):
    if s not in all_states[1]:
        continue
    for n in range(1, 20):
        if s not in all_states[n] or s not in all_states[n+1]:
            continue
        dot = ci_dot(all_states[n][s]['ci'], all_states[n+1][s]['ci'])
        if dot < 0:
            all_states[n+1][s]['ci'] = [
                (o, v, -a) for o, v, a in all_states[n+1][s]['ci']
            ]

with open('all_states_pairphase.pkl', 'wb') as f:
    pickle.dump(all_states, f)
print("Saved all_states_pairphase.pkl")

# Verify: check diagonal signs at 3→4 for all states
from state_overlap import load_swf, build_X, state_overlap
S_MO = load_swf('pairs/frame_3_4/s_wf_3_4_phasefixed')
print("\nS_II(3→4) after pairwise phase fix:")
for s in range(1, 11):
    if s not in all_states[3] or s not in all_states[4]:
        continue
    X3 = build_X(all_states[3][s]['ci'])
    X4 = build_X(all_states[4][s]['ci'])
    S = state_overlap(X3, X4, S_MO)
    print(f"  S_{s}{s}: {S:+.4f}")
