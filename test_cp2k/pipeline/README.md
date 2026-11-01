# CP2K State-Overlap & CI Rotation Pipeline

This directory contains the Python scripts used to compute CIS state overlaps ($S_{IJ} = C_I^T S_{MO} C_J$) and rotate CI vectors between adjacent MD frames using the pseudo-inverse biorthogonal transform.

## Scripts
- `read_molden_mo.py`: Robust sparse Molden reader (pads missing coefficients to 0).
- `build_pair.py`: Builds the Löwdin-normalized MO cross-overlap for two frames.
- `rotate_ci.py`: Fast, norm-preserving CI rotation using `pinv(S_occ) @ X @ pinv(S_virt).T`.
- `state_overlap.py`: Computes CIS overlap matrix elements (vectorized).
- `build_state_overlap_matrix.py`: Loops over trajectories to build $N \times N$ state overlap matrices.
- `phase_fix_*.py`: Handles the arbitrary phase of Davidson eigenvectors for sign-consistent tracking.

## Validation
Tested on 19 adjacent frames of the ML-BP trajectory:
- Norm conservation: `2.5e-5` (machine precision)
- State-1 tracking: `|S_11| = 0.963 - 0.981`
