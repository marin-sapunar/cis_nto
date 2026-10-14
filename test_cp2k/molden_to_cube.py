#!/usr/bin/env python3
"""
Evaluate a single-orbital Molden file on a 3D grid, write a Gaussian cube.

Normalization:
  The cube is normalized so that  ∫ |ψ|² dV = 1  on the grid.
  This puts every cube on the same scale, so a single isovalue (e.g. 0.015
  bohr^-3) produces a comparable isosurface for all orbitals — matching the
  convention used for the publication figures.

Other fixes over the previous version:
  * Origin and spacing are computed once and used both for evaluation and for
    the cube header, so they cannot drift apart.
  * The grid box is built from the atom bounding box + a fixed margin,
    guaranteeing the whole molecule is covered on every axis.
  * Atom block parsing is robust to leading whitespace and to both
    "[Atoms] AU" and "[Atoms] Angs" headers.
  * Coefficients shorter than the AO list are zero-padded instead of
    silently truncating the basis.
"""
import sys
import numpy as np


L_MAP = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}


def cart_terms(l):
    if l == 0:
        return [(0, 0, 0)]
    if l == 1:
        return [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    if l == 2:
        return [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
    if l == 3:
        return [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),
                (1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)]
    if l == 4:
        return [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),
                (1,3,0),(0,3,1),(1,0,3),(0,1,3),(2,2,0),
                (2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)]
    return [(0, 0, 0)]


def read_molden(fname):
    lines = open(fname).readlines()
    atoms = []
    shells = []
    coeffs = []

    i = 0
    while i < len(lines):
        s = lines[i].strip()

        # ---- [Atoms] block --------------------------------------------
        if s.startswith('[Atoms]'):
            ang = 'angs' in s.lower()
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('['):
                p = lines[i].split()
                if len(p) >= 6:
                    Z = int(p[2])
                    x, y, z = float(p[3]), float(p[4]), float(p[5])
                    if ang:
                        x *= 1.8897261246
                        y *= 1.8897261246
                        z *= 1.8897261246
                    atoms.append((Z, x, y, z))
                i += 1
            continue

        # ---- [GTO] block ----------------------------------------------
        if s.startswith('[GTO]'):
            i += 1
            center = 1
            while i < len(lines):
                sl = lines[i].strip()
                if sl.startswith('['):
                    break
                if sl == '':
                    i += 1
                    continue

                p = sl.split()

                if len(p) == 1:
                    try:
                        center = int(p[0])
                        i += 1
                        continue
                    except ValueError:
                        i += 1
                        continue

                if len(p) >= 2 and p[0][0].isalpha():
                    l_char = p[0].lower()
                    n_prim = int(p[1])
                    l = L_MAP.get(l_char, -1)
                    exps, cfs = [], []
                    for _ in range(n_prim):
                        i += 1
                        if i >= len(lines):
                            break
                        pp = lines[i].split()
                        exps.append(float(pp[0]))
                        cfs.append(float(pp[1]))
                    if l >= 0:
                        shells.append((center, l, exps, cfs))
                    i += 1
                    continue

                i += 1
            continue

        # ---- [MO] block -----------------------------------------------
        if s == '[MO]':
            i += 1
            while i < len(lines):
                sl = lines[i].strip()
                if sl.startswith('['):
                    break
                if sl.startswith(('Sym=', 'Ene=', 'Spin=', 'Occup=')):
                    i += 1
                    continue
                p = sl.split()
                if len(p) >= 2:
                    try:
                        coeffs.append(float(p[1]))
                    except ValueError:
                        pass
                i += 1
            continue

        i += 1

    return atoms, shells, coeffs


def main():
    if len(sys.argv) < 3:
        print("Usage: molden_to_cube.py input.molden output.cube [ngrid]")
        sys.exit(1)
    inp = sys.argv[1]
    out = sys.argv[2]
    ngrid = int(sys.argv[3]) if len(sys.argv) > 3 else 60

    atoms, shells, coeffs = read_molden(inp)
    print(f"Atoms: {len(atoms)}  Shells: {len(shells)}  Coeffs: {len(coeffs)}")

    ao = []
    for (center, l, exps, cfs) in shells:
        for term in cart_terms(l):
            ao.append((center, l, term, exps, cfs))
    print(f"Cartesian AOs: {len(ao)}")

    n_ao = len(ao)
    if len(coeffs) < n_ao:
        print(f"WARNING: fewer coeffs ({len(coeffs)}) than AOs ({n_ao}). "
              f"Padding with zeros.")
        coeffs = coeffs + [0.0] * (n_ao - len(coeffs))

    # ---- Grid box: atom bounding box + margin -------------------------
    ax = np.array([a[1] for a in atoms])
    ay = np.array([a[2] for a in atoms])
    az = np.array([a[3] for a in atoms])

    margin = 6.0                              # bohr
    xmin, xmax = ax.min() - margin, ax.max() + margin
    ymin, ymax = ay.min() - margin, ay.max() + margin
    zmin, zmax = az.min() - margin, az.max() + margin

    dx = (xmax - xmin) / (ngrid - 1)
    dy = (ymax - ymin) / (ngrid - 1)
    dz = (zmax - zmin) / (ngrid - 1)

    print(f"Grid origin: ({xmin:.4f}, {ymin:.4f}, {zmin:.4f}) bohr")
    print(f"Grid size  : ({ngrid}, {ngrid}, {ngrid})")
    print(f"Spacing    : dx={dx:.4f}  dy={dy:.4f}  dz={dz:.4f} bohr")

    xg0 = xmin + np.arange(ngrid) * dx
    yg0 = ymin + np.arange(ngrid) * dy
    zg0 = zmin + np.arange(ngrid) * dz

    # ---- Evaluate ------------------------------------------------
    print("Evaluating grid...")
    vals = np.zeros((ngrid, ngrid, ngrid))

    for k in range(n_ao):
        ci, l, (p, q, r), exps, cfs = ao[k]
        c = coeffs[k]
        if c == 0.0:
            continue
        cx, cy, cz = ax[ci-1], ay[ci-1], az[ci-1]
        X, Y, Z = np.meshgrid(xg0 - cx, yg0 - cy, zg0 - cz, indexing='ij')
        r2 = X*X + Y*Y + Z*Z
        rad = np.zeros_like(r2)
        for e, cc in zip(exps, cfs):
            rad += cc * np.exp(-e * r2)
        ang = (X**p) * (Y**q) * (Z**r)
        vals += c * rad * ang
        if k % 100 == 0:
            print(f"  AO {k}/{n_ao}")

    # ---- L2 normalization ----------------------------------------
    # Ensures ∫ |ψ|² dV = 1 on the grid, so a single isovalue produces
    # a comparable isosurface for every cube.
    dV = dx * dy * dz
    norm2 = float(np.sum(vals * vals)) * dV
    if norm2 > 0.0:
        scale = 1.0 / np.sqrt(norm2)
        vals *= scale
        print(f"Normalized: ∫|ψ|² dV before = {norm2:.6e}")
        print(f"  |ψ|_max = {np.max(np.abs(vals)):.6f}   "
              f"(isovalue 0.015 = {100*0.015/np.max(np.abs(vals)):.2f}% of max)")
    else:
        print("WARNING: cube has zero norm — skipping normalization.")

    # ---- Write cube ------------------------------------------------
    print("Writing cube...")
    with open(out, 'w') as f:
        f.write("Cube from Molden\n")
        f.write("MO\n")
        f.write(f"{len(atoms):5d} {xmin:12.6f} {ymin:12.6f} {zmin:12.6f}\n")
        f.write(f"{ngrid:5d} {dx:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write(f"{ngrid:5d} {0.0:12.6f} {dy:12.6f} {0.0:12.6f}\n")
        f.write(f"{ngrid:5d} {0.0:12.6f} {0.0:12.6f} {dz:12.6f}\n")
        for Zi, x, y, z in atoms:
            f.write(f"{Zi:5d} {float(Zi):12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        count = 0
        for ix in range(ngrid):
            for iy in range(ngrid):
                for iz in range(ngrid):
                    f.write(f"{vals[ix,iy,iz]:13.5e}")
                    count += 1
                    if count % 6 == 0:
                        f.write("\n")
        if count % 6 != 0:
            f.write("\n")

    print(f"Wrote {out}")


if __name__ == '__main__':
    main()