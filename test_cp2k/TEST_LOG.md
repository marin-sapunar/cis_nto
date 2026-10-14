# CP2K reader test log

Verified end-to-end. Both test cases pass with the reference values
hard-coded in compare.sh.

## Command

    cd test_cp2k
    bash run_all.sh
    bash compare.sh

    cd mlbp
    bash fetch.sh
    bash run.sh
    bash compare.sh

## Result -- adamantane family

    CP2K test regression check
    ================================================================
    Checking dominant NTO occupations...
      PASS  adamantane   occ=0.9976   (expected 0.9976)
      PASS  1methyladamantane   occ=0.9980   (expected 0.9980)
      PASS  13dimethyladamantane   occ=0.9961   (expected 0.9961)
      PASS  135trimethyladamantane   occ=0.9976   (expected 0.9976)

    Checking cube L2 normalization...
      PASS  s1_occ adamantane   I = 1.000002
      PASS  s1_virt adamantane   I = 1.000002
      PASS  homo adamantane   I = 1.000002
      PASS  lumo adamantane   I = 1.000002
      PASS  s1_occ 1methyladamantane   I = 1.000002
      PASS  s1_virt 1methyladamantane   I = 1.000002
      PASS  homo 1methyladamantane   I = 1.000002
      PASS  lumo 1methyladamantane   I = 1.000002
      PASS  s1_occ 13dimethyladamantane   I = 1.000000
      PASS  s1_virt 13dimethyladamantane   I = 1.000000
      PASS  homo 13dimethyladamantane   I = 1.000000
      PASS  lumo 13dimethyladamantane   I = 1.000000
      PASS  s1_occ 135trimethyladamantane   I = 0.999999
      PASS  s1_virt 135trimethyladamantane   I = 0.999999
      PASS  homo 135trimethyladamantane   I = 0.999999
      PASS  lumo 135trimethyladamantane   I = 0.999999

    ================================================================
      ALL CHECKS PASSED
    ================================================================

## Result -- ML-BP B3LYP

    ML-BP B3LYP regression check
    ================================================================
    Checking dominant NTO occupation...
      PASS  occ=0.9940   (expected 0.9940 +/- 0.005)

    Checking cube L2 normalization...
      PASS  s1_occ   I = 1.000001
      PASS  s1_virt   I = 1.000001
      PASS  homo   I = 1.000001
      PASS  lumo   I = 1.000001

    ================================================================
      ALL CHECKS PASSED
    ================================================================

## What this validates

The reader (cp2k_read.f90 / cp2k_sections.f90) correctly parses:

- CP2K [Atoms] and [GTO] sections in the Molden
- Sparse [MO] coefficient blocks (CP2K omits near-zero AO coefficients)
- TDDFPT Excitation analysis blocks in the CP2K log, with the
  ADDED_MOS virtual-count cap applied when the log references more
  MOs than the Molden contains

## Known upstream issues

Three bugs in existing cis_nto code were found during testing.
They are documented in ISSUES.md and are independent of the CP2K
reader.
