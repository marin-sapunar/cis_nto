# Level-4 Minors (L4M) CISD Overlap Validation

This directory contains the exact CISD wavefunction overlap algorithm using a level-4 minors (L4M) expansion.

## Prerequisites
- gfortran
- OpenBLAS
- Python 3 with pyscf and numpy

## Quick Start
To run the core 18-configuration synthetic validation suite:
make clean
make check_l4m_extended
./check_l4m_extended

## Real-Data Validation (PySCF)
To generate real cross-overlap data for Benzene and Ethylene NAMD steps:
cd scripts && python3 gen_real_data.py && python3 gen_multistate_namd.py && cd ..
make check_real check_multi
./check_real
./check_multi

## Performance Benchmarks
To reproduce the amyloid-scale speedups:
make benchmark_l4m
./benchmark_l4m
