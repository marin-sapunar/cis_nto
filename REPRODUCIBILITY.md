# Reproducibility Guide

This repository accompanies the manuscript:
**"Level-4 Minors for Exact CISD Wavefunction Overlaps: Application to Amyloid Scale Nonadiabatic Dynamics"**

## Quick Start

    git clone https://github.com/hamidzabihi/cis_nto.git
    cd cis_nto/test_cisd
    make check_l4m_extended
    ./check_l4m_extended

Expected output: 18 PASS, 0 FAIL (max error <= 4.3e-18).

## Reproducing Manuscript Tables and Figures

### Table 1: Synthetic Validation (18 configurations)
    cd test_cisd
    make check_l4m_extended
    ./check_l4m_extended

### Table 2: 21-Case Amyloid Validation (77-atom 2B7 geometry)
    cd test_cisd
    python3 gen_amyloid_suite.py
    for f in amyloid_t*_n*.txt; do
        cp "$f" real_data.txt
        ./check_real_tri
    done

Expected: all cases PASS, max error <= 7.3e-12.

### Table 3: Speedup Benchmarks (nstates=16)
    cd test_cisd
    make benchmark_l4m
    ./benchmark_l4m

Expected: 22.46x speedup at nstates=16, n_dets=60.

### Figure 1 & 2: Truncation Convergence and Speedup Plot
    cd test_cisd
    python3 make_plots.py

### Figure 3: ML Schur Surrogate
    cd test_cisd
    python3 prototype_ml_schur.py

### Figure 4: ML Learning Curve
    cd test_cisd
    python3 ml_learning_curve.py

### Python API (15 tests)
    cd python_api
    pytest tests/ -v

## Future Work Prototypes (Discussion Section)

Each prototype directly backs one paragraph in the manuscript.

| File | Claim | Command |
|---|---|---|
| prototype_nto_hybrid.py | 135x singles compression, bounded error | python3 prototype_nto_hybrid.py |
| prototype_mr_cisd_tikhonov_fixed.py | Pseudo-determinant converges to 1.81 | python3 prototype_mr_cisd_tikhonov_fixed.py |
| prototype_krr_schur.py | 8.6e-07 error, O(N^2) memory limit | python3 prototype_krr_schur.py |
| prototype_ann_schur.py | 9.0e-04 error, O(1) inference (23.5 us) | python3 prototype_ann_schur.py |
| prototype_e3nn.py | Zero deviation under 3D rotation | python3 prototype_e3nn.py |

## Dependencies

- Fortran: gfortran, OpenBLAS
- Python: numpy, scipy, scikit-learn, matplotlib, pyscf, torch, e3nn

Install:

    pip install numpy scipy scikit-learn matplotlib pyscf torch e3nn

## Hardware

All benchmarks reported in the manuscript were obtained on a single-core Linux workstation.
