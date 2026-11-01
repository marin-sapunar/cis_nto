import re
import numpy as np

def parse_all_states(log_file):
    with open(log_file) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if "Excitation analysis" in l)
    state_re = re.compile(r'^\s+(\d+)\s+(\d+\.\d+)\s+eV')
    next_re = re.compile(r'^\s+\d+\s+\d+\.\d+\s+eV')
    states = {}
    i = start
    while i < len(lines):
        m = state_re.match(lines[i])
        if m:
            state_num = int(m.group(1))
            energy_ev = float(m.group(2))
            i += 1
            ci = []
            while i < len(lines) and not next_re.match(lines[i]):
                p = lines[i].split()
                if len(p) == 3:
                    try:
                        ci.append((int(p[0]), int(p[1]), float(p[2])))
                    except ValueError:
                        pass
                i += 1
            states[state_num] = {'energy': energy_ev, 'ci': ci}
        else:
            i += 1
    return states

if __name__ == "__main__":
    import pickle
    all_states = {}
    for n in range(1, 21):
        all_states[n] = parse_all_states(f'run_{n}.log')
        n_states = len(all_states[n])
        n_ci = {k: len(v['ci']) for k, v in all_states[n].items()}
        print(f"run_{n}.log: {n_states} states, CI counts = {n_ci}")

    with open('all_states.pkl', 'wb') as f:
        pickle.dump(all_states, f)
    print("\nSaved all_states.pkl")
