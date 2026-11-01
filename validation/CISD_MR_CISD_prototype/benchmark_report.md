# EasyWFO CISD Caching Engine: Performance Benchmark

## Objective
To evaluate the computational performance of the Plasser unique-factor caching algorithm implemented in the C++ EasyWFO CISD bridge for production-scale Non-Adiabatic Molecular Dynamics (NAMD).

## Methodology
We benchmarked the calculation of a \(10 \times 10\) state overlap matrix, where each state consists of 500 CISD determinants. This results in **25,000,000 pairwise determinant overlap evaluations** per time step.

The MO overlap matrix was generated with off-diagonal couplings of up to 0.1. Determinants were generated using a realistic CISD active space (shared core orbitals, with excitations into virtual orbitals).

## Optimization Iterations (Algorithmic)

| Configuration | Data Structure / Compilation | Cache Unique Entries | Execution Time |
| :--- | :--- | :--- | :--- |
| **1. Naive Random Det.** | `std::map` (Red-Black Tree) | 13,645,636 | 21,691 ms |
| **2. Realistic CISD Det.** | `std::unordered_map` (Hash Table) | 23,409 | 647 ms |
| **3. HPC Parallel** | `#pragma omp`, `Matrix4d`, AVX (`-march=native`) | 23,409 | **256 ms** |

## Production Integration Benchmark (Fortran-C++ Bridge)

We ported the HPC optimizations directly into the production `easywfo_cisd_bridge.cpp` used to generate Libra-compatible `st_adi_*.txt` files. We timed the full 50-step trajectory generation (including Fortran I/O):

| Test | Execution Time | Per-Step Time |
| :--- | :--- | :--- |
| **50 NAMD time steps (3 states)** | **0.053 s** | **~1.06 ms** |

## Conclusion
By exploiting the physical structure of CISD wavefunctions (shared core orbitals), utilizing \(O(1)\) hash-map lookups, enforcing static stack memory allocation for the matrices, and parallelizing across 8 CPU threads, the EasyWFO engine achieves **97.6 million overlap evaluations per second** in isolated benchmarks, and **~1.06 milliseconds per full NAMD time step** in production. This performance easily satisfies the requirements for on-the-fly NAMD trajectory simulations.
