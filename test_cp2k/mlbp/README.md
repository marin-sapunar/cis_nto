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
