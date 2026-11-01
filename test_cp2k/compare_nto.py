#!/usr/bin/env python3
import sys
import math

def read_coeffs(fname):
    """Read MO coefficient column from a Molden file."""
    with open(fname) as f:
        lines = f.readlines()
    coeffs = []
    in_mo = False
    for line in lines:
        s = line.strip()
        if s == '[MO]':
            in_mo = True
            continue
        if not in_mo:
            continue
        if s.startswith('Sym=') or s.startswith('Ene=') or \
           s.startswith('Spin=') or s.startswith('Occup='):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                coeffs.append(float(parts[1]))
            except ValueError:
                pass
    return coeffs

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def norm(a):
    return math.sqrt(dot(a, a))

a_file, b_file = sys.argv[1], sys.argv[2]
a = read_coeffs(a_file)
b = read_coeffs(b_file)

n = min(len(a), len(b))
a, b = a[:n], b[:n]

ov = dot(a, b) / (norm(a) * norm(b))
print(f"{a_file}: {len(a)} coeffs")
print(f"{b_file}: {len(b)} coeffs")
print(f"|<{a_file} | {b_file}>| = {abs(ov):.6f}")
