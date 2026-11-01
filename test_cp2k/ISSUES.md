# Upstream bugs found while testing the CP2K reader

Three bugs in existing cis_nto code were found while building the CP2K
tests in this directory. None are caused by the CP2K reader -- the
reader works correctly on all inputs tested. Each should become a
separate GitHub issue.

## Bug 1 -- DGEMM parameter 10 for near-singular CIS vectors

Component: cis_nto.exe (via src/util_chem/l2m/ssblock_lu.c)

Symptom. cis_nto -in cp2k on a TDDFPT calculation whose bright state
is dominated by a single configuration prints many copies of

    ** On entry to DGEMM  parameter number 10 had an illegal value

and produces all-zero NTOs. Other functionals with more mixed
configurations run to completion.

Reproducer. ML-BP / PBEh S1 (99.35% single-config weight):

    cis_nto.exe -in cp2k \
      "PBEh0125_50states/0125_PBE0.molden:PBEh0125_50states/out.log"

Control (passing): ML-BP / B3LYP (mixed S1):

    cis_nto.exe -in cp2k \
      "B3LYP_50states/...molden:B3LYP_50states/out.log"

Input shapes are identical in both cases (nocc=160, nvirt=1344,
nbas=1504). The failure is content-dependent.

Suspected location. src/util_chem/l2m/ssblock_lu.c:88:

    cblas_dgemm(..., nv1, no/2, nv2, zero, wf1, nv1, wf2, nv2,
                neg_one, &csc[no+nv1 + no], 2*ld_csc);

The 2*ld_csc stride starting at &csc[no+nv1+no] reads or writes
outside the csc allocation for certain coefficient patterns.

## Bug 2 -- convert.exe -out molden_cart sparse-write corruption

Component: convert.exe (via molden_read.f90 / molden_write.f90)

Symptom. Converting a CP2K Molden with sparse coefficient blocks
to Cartesian form loses most coefficients for some MOs.

Reproducer. ML-BP / PBEh MO 160. In std_PBEh.molden (output of
convert.exe -in cp2k -out molden_cart), the MO has 1 nonzero
coefficient out of 896. The source Molden
(PBEh0125_50states/0125_PBE0.molden) has 798 nonzero coefficients
for that MO. Deterministic for a given input.

Consequence. Cubes generated from the standardized Molden for that
MO are physically wrong -- density concentrated on a single atom.

Workaround. Extract the MO directly from the source Molden with a
dense extractor, bypassing convert.exe.

## Bug 3 -- NTO Molden writer sparse-coefficient corruption

Component: cis_nto.exe (NTO Molden writer path)

Symptom. For some state/NTO combinations, the block written by
cis_nto to nto.molden has near-zero coefficient amplitudes even
though the Occup= header reports the correct (high) occupation.

Reproducer. ML-BP / B3LYP. In nto_B3LYP.molden, block
Sym= 2001a (Occup= 0.994035) has max|coeff| = 0.163 -- an order
of magnitude smaller than a typical MO coefficient. The companion
block Sym= 1001a (Occup= 0.989472) is correct.

The block header is fine; only the coefficient values are wrong.

Consequence. Cubes built from that block are degenerate (very
small raw grid norm). After grid normalization they integrate to 1,
so a naive norm check does not catch the corruption.
