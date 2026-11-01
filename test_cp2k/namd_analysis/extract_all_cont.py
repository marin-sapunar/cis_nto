import re, pickle
def parse_all_states(log_file):
    with open(log_file) as f: lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if "Excitation analysis" in l)
    state_re = re.compile(r'^\s+(\d+)\s+(\d+\.\d+)\s+eV')
    next_re = re.compile(r'^\s+\d+\s+\d+\.\d+\s+eV')
    states = {}
    i = start
    while i < len(lines):
        m = state_re.match(lines[i])
        if m:
            s = int(m.group(1)); e = float(m.group(2)); i += 1; ci = []
            while i < len(lines) and not next_re.match(lines[i]):
                p = lines[i].split()
                if len(p) == 3:
                    try: ci.append((int(p[0]), int(p[1]), float(p[2])))
                    except ValueError: pass
                i += 1
            states[s] = {'energy': e, 'ci': ci}
        else: i += 1
    return states

for n in range(1, 21):
    try:
        states = parse_all_states(f'run_cont_{n}.log')
        with open(f'states_cont_{n}.pkl', 'wb') as f: pickle.dump(states, f)
    except Exception as ex: print(f"Failed on {n}: {ex}")
print("All 20 states extracted.")
