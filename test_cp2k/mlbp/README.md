# ML-BP test case (fetched from upstream)

This is a second, larger test for the CP2K TDDFPT reader, using the
ML-BP (multilayer black phosphorus) dataset published in:

    https://github.com/AkimovLab/Project_MLBP_excitonic

The complete dataset contains seven exchange-correlation functionals
(PBEh, PBE0, B3LYP, BHandHLYP, CAM-B3LYP, wB97X, wB97XD), each with 50
TDDFPT excited states. To keep the `cis_nto` repository small and avoid
redistributing third-party data, this test fetches only one functional
(B3LYP) from the upstream repository at run time.

## What it exercises

The ML-BP system is a 64-atom periodic supercell (832 SCF molecular
orbitals, 160 occupied). It exercises several code paths that the
adamantane test does not:

- **Sparse Molden coefficient blocks** — CP2K omits near-zero AO
  coefficients. The reader places each coefficient at its true AO
  index instead of assuming a dense sequence.
- **`ADDED_MOS` virtual mismatch** — CP2K ≤ 2024.1 does not write the
  `&MO_MOLDEN ADDED_MOS` virtuals to the Molden. The reader caps
  `added_mos` at the Molden's actual MO count and issues a warning.
- **Large system** — 64 atoms, 896 Cartesian AOs, 1504 total MOs in
  the CIS space.

## Running

    bash fetch.sh    # download the B3LYP Molden + log from AkimovLab
    bash run.sh      # run the CP2K pipeline
    bash compare.sh  # verify the output

`fetch.sh` uses `curl` to download two files:

    b3lyp_opted-MOS-1_0.molden   (~16 MB)
    out.log                      (~1 MB)

Both are copied from the `B3LYP_50states/` directory of the upstream
repository at `main`.

## Expected output

For B3LYP, the dominant S1 NTO pair should have occupation ~0.99
(consistent with the ~98% single-configuration character of the S1
state). Every generated cube should satisfy ∫|ψ|² dV = 1.000000.

## Upstream citation

When publishing results obtained with this test, cite the ML-BP
dataset:

    [Citation to be filled in by the user from the AkimovLab repo's
     README or associated publication.]