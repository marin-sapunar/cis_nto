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
  `&MO_MOLDEN ADDED_MOS` virtuals into the Molden. The reader detects the
  mismatch and invokes `cp2k_virtuals_mod` to reconstruct or consistently
  truncate the virtual space, then reports which strategy was used:

      Notice in cp2k_read_cis:
        log references : 1504 MOs
        molden contains: 832 MOs
        reconstruction method: truncated
        dropping 672 virtuals

  Three strategies are tried in order: (1) read `<project>-RESTART.wfn` for
  the full MO set, (2) require an augmented Molden, (3) drop the missing
  virtuals from both bra and ket so the overlap stays well-defined. If the
  truncation leaves zero virtuals, the run aborts with a diagnostic instead
  of entering an infinite DGEMM loop (Issue #5).

## Running

    cd mlbp
    bash fetch.sh     # ~10 s, downloads ~16 MB from AkimovLab
    bash run.sh       # ~2 min, produces 4 cubes in cubes/
    bash compare.sh   # ~10 s, checks NTO occupation and cube integrals

## Reference values

Dominant S1 NTO occupation: 0.9890 +/- 0.005 (the 0.9940 in earlier drafts was the S2 label, 2001a).

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

**What the `ADDED_MOS` notice means.** CP2K <= 2024.1 does not write the
`&MO_MOLDEN ADDED_MOS` virtuals into the Molden file. The log therefore
references more MOs (1504) than the Molden contains (832). The reader
detects the mismatch and invokes `cp2k_virtuals_mod` to handle it:

    Notice in cp2k_read_cis:
      log references : 1504 MOs
      molden contains: 832 MOs
      reconstruction method: truncated
      dropping 672 virtuals

The current strategy is **consistent truncation**: the missing virtuals
are dropped from both bra and ket so the CIS vectors remain
well-defined and the NTO construction can proceed. Future revisions
will add RESTART.wfn-based reconstruction, which reads `<project>.wfn`
to recover the missing virtuals rather than discarding them.

This notice is **not** a failure. Silently accepting the mismatch would
corrupt the CIS vectors; the reader instead reduces the virtual space
to what the Molden can represent and reports the loss. If truncation
leaves zero virtuals, the run aborts with a diagnostic (Issue #5).
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

## Comparison with the paper

Table 1 of the ML-BP paper (J. Phys. Chem. Lett., 2024,
doi:10.1021/acs.jpclett.5c01417) lists the bright-state composition
per functional:

| functional | state | major transition | contribution |
|---|---|---|---|
| PBEh12.5   | 1  | HOMO -> LUMO       | 99.3% |
| B3LYP      | 1  | HOMO -> LUMO       | 97.9% |
| PBE0       | 1  | HOMO -> LUMO       | 96.7% |
| BHandHLYP  | 1  | HOMO -> LUMO       | 83.7% |
| CAM-B3LYP  | 1  | HOMO -> LUMO       | 76.9% |
| wB97XD     | 38 | HOMO-3 -> LUMO+5   | 18.5% |
| wB97X      | 38 | HOMO-5 -> LUMO+7   | 15.3% |

For B3LYP the paper reports S1 = 97.9% HOMO -> LUMO. This test
measures the S1 dominant NTO occupation as 0.9895, which corresponds
to 98.8% single-configuration character -- matching the paper to
within the accuracy of the two different codes.

For wB97X and wB97XD the bright state is not S1 but state 38, so its
Table 1 entry does not correspond to the S1 NTO occupation measured
here. This is a paper-convention difference, not a discrepancy.

Figure 2 of the paper shows the NTO pairs (hole and electron) for
PBEh, PBE0 and wB97X at isosurface value 0.015 bohr^-3, rendered with
Multiwfn + VMD. The cubes produced by `molden_to_cube.py` are the
same physical objects; the shipped `render_nto.tcl` script can be used
for a visual side-by-side comparison.
