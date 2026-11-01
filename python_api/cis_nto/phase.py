"""Phase correction and state tracking for overlap matrices.

Implements the Local Diabatization (LD) phase-tracking approach:
1. Hungarian permutation to align states across time steps
2. Sign flip so diagonal overlaps are positive real

Reference: Shakiba & Akimov (2023), Generalization of the local
diabatization approach for propagating electronic DOF in NAMD.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment


def correct_phase_matrix(overlap, return_swaps=False):
    """Correct sign flips and state permutations in overlap matrices.

    Parameters
    ----------
    overlap : np.ndarray (nsteps, nstates, nstates) complex
        Overlap matrices S_IJ(t) = <Psi_I(t) | Psi_J(t+dt)>
    return_swaps : bool
        If True, also return per-step swap and flip counts.

    Returns
    -------
    corrected : np.ndarray (nsteps, nstates, nstates)
    swaps : list[int]  (if return_swaps)
    flips : list[int]  (if return_swaps)
    """
    overlap = np.asarray(overlap, dtype=complex)
    if overlap.ndim != 3:
        raise ValueError(f"Expected 3D overlap (nsteps,nstates,nstates), got {overlap.shape}")
    nsteps, nstates, ncols = overlap.shape
    if nstates != ncols:
        raise ValueError(f"Last two dims must be equal, got {nstates}x{ncols}")

    corrected = overlap.copy()
    swaps, flips = [], []

    for t in range(nsteps):
        S = corrected[t]
        # Hungarian algorithm: minimize negative |S_ij| = maximize |S_ij|
        row_ind, col_ind = linear_sum_assignment(-np.abs(S))
        perm = np.zeros(nstates, dtype=int)
        perm[row_ind] = col_ind
        S_perm = S[:, perm]
        # Sign flip: diagonal must be positive real
        diag_real = np.real(np.diag(S_perm))
        signs = np.where(diag_real < 0, -1.0, 1.0)
        S_corrected = S_perm * signs[np.newaxis, :]
        corrected[t] = S_corrected
        swaps.append(int(np.sum(perm != np.arange(nstates))))
        flips.append(int(np.sum(signs < 0)))

    if return_swaps:
        return corrected, swaps, flips
    return corrected


def verify_corrected(overlap, threshold=0.5):
    """Verify a corrected overlap has all-positive diagonal above threshold.

    Returns
    -------
    ok : bool
    n_bad : int
        Number of diagonal elements below threshold
    """
    overlap = np.asarray(overlap, dtype=complex)
    diag = np.real(np.diagonal(overlap, axis1=1, axis2=2))
    n_bad = int(np.sum(diag < threshold))
    return n_bad == 0, n_bad
