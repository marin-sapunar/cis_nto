# NAMD Analysis Scripts

These scripts demonstrate that the CP2K overlap pipeline is ready for Non-Adiabatic Molecular Dynamics (FSSH).

## Scripts
- `extract_all_cont.py`: Parses the `Excitation analysis` block from CP2K logs into a pickle file of CI vectors.
- `calculate_nacs.py`: Computes the Time-Derivative Non-Adiabatic Couplings (NACs) and hopping probabilities using the pseudo-inverse overlap matrix.
- `test_overlaps_cont.py`: A quick test to compute \(S_{11}\) tracking between two adjacent frames.
- `plot_pes.py`: Plots the Potential Energy Surfaces (PES) for all 15 states across the trajectory.
- `run_fssh.py`: A toy FSSH script demonstrating surface hopping using the calculated probabilities.

## Prerequisites
These scripts rely on `state_overlap.py` and `read_molden_mo.py` located in `test_cp2k/pipeline/`.
