# HPC-Optimized CISD/MR-CISD Fortran-C++ Prototype

This directory contains the production-ready prototype for calculating CISD/MR-CISD wavefunction overlaps. It uses the Plasser unique-factor caching algorithm implemented in C++ and interfaced with Fortran.

## Performance Highlights
- **Algorithmic Speedup:** 84.7x (21.6s down to 256ms for 25M evaluations) via `std::unordered_map` and `Matrix4d` stack allocation.
- **Production Speed:** ~1.06 ms per full NAMD time step (generates a 50-step trajectory in 0.051s).
- **Robustness:** Base-10000 hash keys support active spaces up to 9,999 orbitals with zero collisions.
- **Parallelization:** Utilizes OpenMP to distribute the state overlap matrix calculation across all available CPU cores.

## Reproducibility Steps

### 1. Generate the Libra Trajectory (Production Test)
    make clean
    make
    ./test_fortran_trajectory

### 2. Run the HPC Benchmark
    g++ -O3 -march=native -fopenmp -I /usr/include/eigen3 bench_cisd_hpc.cpp -o bench_cisd_hpc
    ./bench_cisd_hpc

### 3. Run the Large Active Space Stress Test
    g++ -O3 -I /usr/include/eigen3 stress_test_large_mo.cpp -o stress_test_large_mo
    ./stress_test_large_mo

## Documentation
See `benchmark_report.md` for the full optimization journey (21.6s → 647ms → 256ms).
