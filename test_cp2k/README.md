# CP2K reader test cases

Reproducible tests for the CP2K TDDFPT reader added to `cis_nto`:

- `cp2k_read.f90` + `cp2k_sections.f90` — new readers
- `molden_read.f90` + `molden_sections.f90` — Molden parsing (with
  sparse-coefficient fix)
- `read_all.f90` — dispatcher gains `case('cp2k')`

There are two test cases in this directory:

1. **Adamantane family** (in this directory) — four small molecules,
   dense Molden output. Tests basic correctness.
2. **ML-BP B3LYP** (in `mlbp/`) — a 64-atom periodic supercell with
   sparse Molden output, exercising the `ADDED_MOS` mismatch and larger
   systems.

See `mlbp/README.md` for the ML-BP test.

---

## Test case 1: Adamantane family

Four adamantane-family molecules of increasing size. Each has a Molden
file (geometry, basis, MOs) and a CP2K log (TDDFPT excitations).

### System

Adamantane family, C10H16 and methylated derivatives:
- `adamantane` — 26 atoms, 68 MOs (38 occupied + 30 virtual)
- `1_methyladamantane` — 29 atoms, 72 MOs
- `1_3_dimethyladamantane` — 32 atoms, 76 MOs
- `1_3_5_trimethyladamantane` — 35 atoms, 80 MOs

GTH pseudopotentials, DZVP-MOLOPT-GTH basis. ECP valence-only
calculations.

### Path convention

The CP2K reader accepts a compound path string:

    <molden_path>:<log_path>

Both files are required because CP2K writes MO coefficients to a Molden
file and TDDFPT excitations to the log. The two halves are split on the
first `:`.

### Files in this directory

Sources (inputs):
- `adamantane.molden`, `adamantane.log`
- `1_methyladamantane.molden`, `1_methyladamantane.log`
- `1_3_dimethyladamantane.molden`, `1_3_dimethyladamantane.log`
- `1_3_5_trimethyladamantane.molden`, `1_3_5_trimethyladamantane.log`

Scripts:
- `run_all.sh` — driver: converts moldens, runs `cis_nto`, extracts
  NTO pairs, generates HOMO/LUMO and S1 cube files
- `compare.sh` — regression check: verifies NTO occupations and cube
  L2 norms against reference values
- `extract_nto.py` — extracts a single Sym-tagged MO block from a
  Molden file (exact Sym match, handles both Sym-delimited and
  Ene-delimited dialects)
- `molden_to_cube.py` — evaluates a single-MO Molden file on a 3D grid
  and writes a Gaussian cube (with L2 normalization)
- `compare_nto.py` — computes overlap integrals between two Molden
  MO coefficient vectors
- `make_grid.py` — helper for grid-based analysis
- `render_all.tcl`, `render_nto.tcl` — VMD scripts for isosurface
  rendering

### Running

From this directory:

    bash run_all.sh
    bash compare.sh

`run_all.sh` for each molecule:

1. Converts the CP2K Molden to a standardized form via `convert.exe`
2. Runs `cis_nto.exe -in cp2k` on the raw CP2K input, producing
   `nto_<tag>.molden`
3. Extracts the dominant S1 NTO pair into `s1_occ_<tag>.molden` and
   `s1_virt_<tag>.molden`
4. Extracts the S0 HOMO/LUMO into `homo_<tag>.molden` and
   `lumo_<tag>.molden`
5. Generates cube files for all four orbitals under `cubes/`

`compare.sh` then checks:

- The dominant S1 NTO occupation for each molecule is within 1e-3 of
  the expected value (0.996–0.998).
- Each generated cube satisfies ∫|ψ|² dV = 1 to within 1e-3.

### Expected output

For each molecule, `cis_nto.exe` reports the number of excited states,
orbitals, and occupied/virtual partitions that match the CP2K log.
For adamantane:

    Number of excited states:           30
    Number of orbitals:                 98
    Number of occupied orbitals:        38
    Number of virtual orbitals:         60

(The virtual count includes CP2K's `ADDED_MOS` virtuals, which are not
written to the Molden file by default in CP2K ≤ 2024.1. The reader
detects this mismatch and caps `added_mos` at the Molden's actual MO
count, printing a warning.)

The dominant NTO pair has occupation ~0.9976 for every molecule in this
set, consistent with the S1 state being dominated by a single
HOMO→LUMO transition.

### Correctness indicators

Two independent checks confirm the interface is correct:

1. **NTO occupations**: the dominant pair has occupation > 0.99,
   the rest have occupation < 0.01. A clean bimodal distribution.
   If the AO ordering is misread, occupations become spread across
   0.1–0.9 with no clean separation.

2. **Cube normalization**: every generated cube satisfies
   ∫|ψ|² dV = 1.000000 to numerical precision. If the AO-to-coefficient
   mapping is wrong, this integral diverges from 1.

---

## Test case 2: ML-BP B3LYP (in `mlbp/`)

A second, larger test uses the ML-BP dataset from
[AkimovLab/Project_MLBP_excitonic](https://github.com/AkimovLab/Project_MLBP_excitonic).

The ML-BP input files are large (a single functional's Molden file is
~16 MB) and come from an external repository without an explicit
redistribution license. To keep this repository small and avoid
redistributing third-party data, the ML-BP test downloads the B3LYP
input at run time:

    cd mlbp
    bash fetch.sh     # downloads b3lyp_opted-MOS-1_0.molden + out.log
    bash run.sh       # runs the pipeline
    bash compare.sh   # verifies the output

The ML-BP test exercises:

- **Sparse Molden coefficient blocks** — CP2K omits near-zero AO
  coefficients in its default output; the reader places each
  coefficient at its true AO index instead of assuming a dense
  sequence.
- **`ADDED_MOS` virtual mismatch** — CP2K ≤ 2024.1 does not write the
  `&MO_MOLDEN ADDED_MOS` virtuals to the Molden; the reader caps
  `added_mos` at the Molden's actual MO count and issues a warning.
- **Large system** — 64 atoms, 832 SCF molecular orbitals, 1504 total
  MOs in the CIS space (672 added virtuals).

The complete seven-functional dataset (PBEh, PBE0, B3LYP, BHandHLYP,
CAM-B3LYP, wB97X, wB97XD) is available at the upstream repository.

---

## Known limitations

- Unrestricted (spin-polarized) TDDFPT is not implemented. Only
  restricted closed-shell calculations are supported.
- CP2K ≤ 2024.1 does not write `&MO_MOLDEN ADDED_MOS` virtuals to the
  Molden file. When the TDDFPT section requests more virtuals than the
  Molden provides, the reader caps `added_mos` at the Molden's MO count
  and issues a warning.

## Files added to `cis_nto`

Two new files in `src/util_chem/io/`:

- `cp2k_read.f90` — public entry points (`cp2k_read_geom`,
  `cp2k_read_basis`, `cp2k_read_mo`, `cp2k_read_cis`)
- `cp2k_sections.f90` — low-level parsers for Molden and log sections

Three modified files:

- `read_all.f90` — added `case('cp2k')` branches and a bounds-safe
  format-string selector
- `molden_sections.f90` — sparse coefficient handling (CP2K omits
  near-zero coefficients from the Molden output)
- `cp2k_read.f90` — `added_mos` cap and two-column spin masks