# CP2K reader test cases

Reproducible tests for the CP2K TDDFPT reader added to `cis_nto`.

## What this adds

New files in `src/util_chem/io/`:

- `cp2k_read.f90`      - public entry points (`cp2k_read_geom`,
                         `cp2k_read_basis`, `cp2k_read_mo`, `cp2k_read_cis`)
- `cp2k_sections.f90`  - low-level parsers for CP2K Molden and log

Modified:

- `read_all.f90`        - `case('cp2k')` branches, bounds-safe format selector
- `molden_sections.f90` - sparse coefficient handling

## Test cases

1. **Adamantane family** (this directory) - four small molecules with
   dense Molden output.
2. **ML-BP B3LYP** (`mlbp/`) - a 64-atom periodic supercell with
   sparse Molden output.

---

## Test 1 - Adamantane family

### System

| tag | atoms | MOs | nocc | nvirt |
|---|---|---|---|---|
| `adamantane` | 26 | 68 | 38 | 30 |
| `1methyladamantane` | 29 | 72 | 38 | 34 |
| `13dimethyladamantane` | 32 | 76 | 38 | 38 |
| `135trimethyladamantane` | 35 | 80 | 38 | 42 |

GTH pseudopotentials, DZVP-MOLOPT-GTH basis, ECP valence-only.

### Path convention

The CP2K reader accepts a compound path string:

    <molden_path>:<log_path>

Both files are required because CP2K writes MO coefficients to a Molden
file and TDDFPT excitations to the log. The two halves are split on the
first `:`.

### Running

    cd test_cp2k
    bash run_all.sh     # ~40 s, produces 16 cubes in cubes/
    bash compare.sh     # ~5 s, checks NTO occupations and cube integrals

### Pipeline

`run_all.sh` per molecule:

1. `convert.exe -in cp2k -out molden_cart` -> `std_<tag>.molden`
2. `cis_nto.exe -in cp2k` on the raw CP2K pair -> `nto_<tag>.molden`
3. Extracts the S1 NTO pair (dominant occupied + dominant virtual)
4. Extracts HOMO/LUMO by block index (38, 39) from `std_<tag>.molden`
5. Generates cubes for all four orbitals

### Reference values

Dominant S1 NTO occupation:

| molecule | occ |
|---|---|
| adamantane | 0.9976 |
| 1methyladamantane | 0.9980 |
| 13dimethyladamantane | 0.9961 |
| 135trimethyladamantane | 0.9976 |

Each cube satisfies the L2 integral = 1.000000 on the grid, with at
least 2000 non-negligible cells.

### Correctness indicators

1. NTO occupations - clean bimodal, dominant > 0.99.
2. Cube integrals - L2 norm = 1 to 1e-3.
3. Nonzero cell count - floor of 2000 catches near-empty MO blocks
   that would still renormalize to 1.

---

## Test 2 - ML-BP B3LYP

Uses the ML-BP dataset from
[AkimovLab/Project_MLBP_excitonic](https://github.com/AkimovLab/Project_MLBP_excitonic).
Only B3LYP is fetched at run time, to keep the repository small and
avoid redistributing third-party data.

    cd mlbp
    bash fetch.sh
    bash run.sh
    bash compare.sh

Exercises sparse Molden coefficient blocks, the `ADDED_MOS` virtual
mismatch, and a large system (64 atoms, 832 SCF MOs, 1504 total CIS
MOs, 896 Cartesian AOs).

---

## Supporting files

- `extract_nto.py` - extracts a single MO block by index or Sym tag
  (refuses ambiguous Sym matches)
- `molden_to_cube.py` - evaluates a single-MO Molden on a 3D grid and
  writes a Gaussian cube (with L2 normalization)
- `validate_against_paper.py` - standalone tool: parses the CP2K log's
  `Excitation analysis` block and compares a state's leading transition
  amplitude against a user-supplied reference value. Not called by
  `compare.sh`.
- `compare_nto.py` - overlap between two Molden coefficient vectors
- `make_grid.py` - grid helper
- `render_all.tcl`, `render_nto.tcl` - VMD rendering
- `TEST_LOG.md` - output of a passing run
- `RUN_TESTS.md` - quick-start guide
- `ISSUES.md` - upstream bugs found during testing

## Known limitations

- Unrestricted (spin-polarized) TDDFPT is not supported.
- CP2K <= 2024.1 does not write `&MO_MOLDEN ADDED_MOS` virtuals into
  the Molden. The reader caps `added_mos` at the Molden's MO count
  and warns.
- `convert.exe -out molden_cart` and the NTO Molden writer both have
  sparse-coefficient issues for some CP2K inputs. See `ISSUES.md` for
  reproduction cases. The tests here avoid the affected MOs by using
  explicit block indices.
