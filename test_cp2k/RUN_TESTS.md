# Running the CP2K reader tests

Two test cases. Both finish in about a minute on a normal workstation.

## Requirements

- cis_nto built with the CP2K reader (default CIS_NTO_BIN=/home/hamid/cis_nto/build/bin)
- python3 with numpy
- curl for the ML-BP fetch
- bash, awk

Override the binary directory with:

    export CIS_NTO_BIN=/path/to/cis_nto/build/bin

## Test 1 -- Adamantane family

    cd test_cp2k
    bash run_all.sh
    bash compare.sh

run_all.sh produces 16 cubes in cubes/. compare.sh checks NTO
occupations and cube integrals.

Expected last line: ALL CHECKS PASSED.

Runtime: ~40 s.

## Test 2 -- ML-BP B3LYP

    cd test_cp2k/mlbp
    bash fetch.sh
    bash run.sh
    bash compare.sh

fetch.sh downloads ~16 MB from the AkimovLab repository. run.sh
produces 4 cubes. compare.sh checks the NTO occupation and cube
integrals.

Expected last line: ALL CHECKS PASSED.

Runtime: ~2 min (mostly grid evaluation).

## If a check fails

compare.sh prints FAIL and the failing value.

- nto.molden missing -> cis_nto failed. Check nto_<tag>.log.
- Cube integral far from 1.0 -> the Molden was misread. Confirm
  cis_nto was built from the same source tree as the running binary.
- nonzero below the floor -> the extracted MO block was almost
  empty. This can happen if the source Molden has sparse-write
  corruption; see ISSUES.md.

See TEST_LOG.md for the exact output of a passing run.
