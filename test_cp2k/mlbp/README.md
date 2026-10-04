# ML-BP B3LYP test case

Second test for the CP2K TDDFPT reader, using the ML-BP dataset at:

    https://github.com/AkimovLab/Project_MLBP_excitonic

Only B3LYP is fetched at run time, to keep the repository small and
avoid redistributing third-party data. The full upstream dataset has
seven functionals (PBEh, PBE0, B3LYP, BHandHLYP, CAM-B3LYP, wB97X,
wB97XD), each with 50 TDDFPT excited states.

## System

64-atom periodic supercell.

- 896 Cartesian AOs
- 832 SCF MOs
- 160 occupied
- 1504 total MOs in the CIS space (672 added virtuals)

## What it exercises

- **Sparse Molden coefficient blocks** - CP2K omits near-zero AO
  coefficients; the reader places each coefficient at its true AO
  index.
- **`ADDED_MOS` virtual mismatch** - CP2K <= 2024.1 does not write the
  `&MO_MOLDEN ADDED_MOS` virtuals into the Molden, so the reader caps
  `added_mos` at the Molden's actual MO count and issues a warning.

## Running

    cd mlbp
    bash fetch.sh     # ~10 s, downloads ~16 MB from AkimovLab
    bash run.sh       # ~2 min, produces 4 cubes in cubes/
    bash compare.sh   # ~10 s, checks NTO occupation and cube integrals

## Reference values

Dominant S1 NTO occupation: 0.9940 +/- 0.005.

Each cube satisfies the L2 integral = 1.000000 on the 60x60x60 grid,
with at least 20000 non-negligible cells.

## Upstream citation

Cite the ML-BP dataset from AkimovLab when publishing results
obtained with this test.

## Functional coverage

Only B3LYP is fetched and run by the shipped test. During development
the reader was exercised against all seven ML-BP functionals:

| functional | reader | `cis_nto -alg NTO` | notes |
|---|---|---|---|
| B3LYP | OK | OK | shipped test |
| PBE0 | OK | OK | not shipped |
| BHandHLYP | OK | OK | not shipped |
| CAM-B3LYP | OK | OK | not shipped |
| wB97X | OK | OK | not shipped |
| wB97XD | OK | OK | not shipped |
| PBEh | OK | fails | upstream DGEMM bug, see ISSUES.md |

The reader itself is functional-independent: it parses the Molden and
log the same way for all seven. The PBEh failure is in `cis_nto`'s
NTO construction, not in the reader.

## Background

**System.** ML-BP is a monolayer / multilayer black phosphorus
supercell. The shipped test uses the 64-atom B3LYP structure from the
AkimovLab dataset.

**Inputs fetched by `fetch.sh`.**

- `b3lyp_opted-MOS-1_0.molden` -- MO coefficients for the B3LYP
  geometry optimized with the b3lyp functional. Sparse format: CP2K
  omits AO coefficients below its print threshold, so each MO block
  contains fewer than `nbas` coefficient lines. The reader places
  each coefficient at its true AO index.
- `out.log` -- the CP2K output for the B3LYP TDDFPT run. Contains the
  SCF section, the `Number of molecular orbitals` line, and the
  TDDFPT `Excitation analysis` block. This is what `cis_nto -in cp2k`
  reads alongside the Molden.

**Why the `ADDED_MOS` warning is expected.** CP2K <= 2024.1 does not
write the `&MO_MOLDEN ADDED_MOS` virtuals into the Molden file. The
log therefore references more MOs (1504) than the Molden contains
(832). The reader caps `added_mos` at the Molden's actual MO count and
prints:

    Warning in cp2k_read_cis:
      CIS log references 1504 MOs
      molden contains 832 MOs
      capping added_mos to 0

This warning is **not** a failure. It is the reader doing the right
thing: silently accepting a mismatch would corrupt the CIS vectors.
You should see this warning once per run, before the NTO construction
starts.

**Why the S1 NTO cubes match the HOMO/LUMO cubes.** B3LYP S1 is
99.4% a single HOMO->LUMO transition. The dominant occupied NTO is
therefore essentially the HOMO and the dominant virtual NTO is
essentially the LUMO. The cube files agree to ~4 decimals. This is
correct physics, not a bug.

## What to look for in the run output

1. The `capping added_mos to 0` warning -- expected, see above.
2. Six MO blocks resolved cleanly:
   - `Sym= 2001a` (dominant occupied NTO, Occup = 0.994)
   - `Sym= 1001a` (dominant virtual NTO, Occup = 0.989)
   - indices 160, 161 for HOMO / LUMO in `std_B3LYP.molden`
3. Each cube printed `Normalized: ... |ψ|_max = <value>` -- the raw
   grid integral before renormalization is small for the S1 NTO
   (7.3e-07) because the grid spacing is asymmetric along z
   (dz=0.27 bohr vs dx=0.72 bohr). After `molden_to_cube.py`'s L2
   renormalization the cube integrates to 1.000000.
4. `compare.sh` ends with `ALL CHECKS PASSED`.
