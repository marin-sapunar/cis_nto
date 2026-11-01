import numpy as np

def read_molden_mo(path):
    """
    Read MO coefficients from a molden file, padding sparse blocks to n_ao.

    CP2K with NDIGITS 9 omits coefficients below 1e-9, so some MOs have
    fewer than n_ao entries. This reader uses the index column to place
    each coefficient at the correct AO position, filling missing entries
    with 0.
    """
    with open(path) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == '[MO]')
    coefs = []
    n_ao = None
    i = start + 1
    while i < len(lines):
        if lines[i].lstrip().startswith('Ene='):
            i += 3
            entries = {}
            max_idx = 0
            while i < len(lines) and not lines[i].lstrip().startswith('Ene=') \
                             and not lines[i].lstrip().startswith('['):
                p = lines[i].split()
                if len(p) >= 2:
                    try:
                        idx = int(p[0])
                        entries[idx - 1] = float(p[1])
                        if idx > max_idx:
                            max_idx = idx
                    except (ValueError, IndexError):
                        pass
                i += 1
            if n_ao is None:
                n_ao = max_idx
            vec = np.zeros(n_ao)
            for k, v in entries.items():
                if 0 <= k < n_ao:
                    vec[k] = v
            coefs.append(vec)
        else:
            i += 1
    return np.array(coefs).T
