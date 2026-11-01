"""Brute-force CISD overlap calculator (reference implementation).

Purpose: ground truth for verifying the fast CISD overlap algorithms.

Approach: explicit determinant-by-determinant summation. O(N_dets^2).
Slow but provably correct.

Reference: Slater-Condon rules for determinant overlaps.
"""
import numpy as np
from itertools import combinations
from typing import List, Tuple, Dict


Det = Tuple[Tuple[int, ...], Tuple[int, ...]]  # (alpha_occ, beta_occ)


def make_determinant(alpha_occ, beta_occ) -> Det:
    """Canonicalize a determinant: sorted occupied orbital indices."""
    return (tuple(sorted(alpha_occ)), tuple(sorted(beta_occ)))


def gen_cisd_determinants(n_occ, n_virt, n_electrons_alpha=1, n_electrons_beta=1):
    """Generate all CISD determinants for a small reference.

    For a minimal reference with n_electrons_alpha occupied alpha and
    n_electrons_beta occupied beta, generate:
    - 1 reference
    - singles (1 alpha OR 1 beta promoted)
    - doubles (1 alpha AND 1 beta promoted)

    Args:
        n_occ: number of occupied orbitals (reference)
        n_virt: number of virtual orbitals
        n_electrons_alpha, n_electrons_beta: reference electron counts

    Returns: (dets, labels) where labels[i] is a string tag.
    """
    # Reference: fully occupied (0..n_occ-1) in alpha and beta
    ref_alpha = tuple(range(n_occ))
    ref_beta = tuple(range(n_occ))
    dets = [(ref_alpha, ref_beta)]
    labels = ["ref"]

    # Singles: promote alpha_i -> virt_a
    for i in range(n_occ):
        for a in range(n_occ, n_occ + n_virt):
            new_alpha = tuple(sorted([o for o in ref_alpha if o != i] + [a]))
            dets.append((new_alpha, ref_beta))
            labels.append(f"alpha_{i}->{a}")

    # Singles: promote beta_i -> virt_a
    for i in range(n_occ):
        for a in range(n_occ, n_occ + n_virt):
            new_beta = tuple(sorted([o for o in ref_beta if o != i] + [a]))
            dets.append((ref_alpha, new_beta))
            labels.append(f"beta_{i}->{a}")

    # Doubles: alpha_i -> virt_a AND beta_j -> virt_b
    for i in range(n_occ):
        for j in range(n_occ):
            for a in range(n_occ, n_occ + n_virt):
                for b in range(n_occ, n_occ + n_virt):
                    new_alpha = tuple(sorted([o for o in ref_alpha if o != i] + [a]))
                    new_beta = tuple(sorted([o for o in ref_beta if o != j] + [b]))
                    dets.append((new_alpha, new_beta))
                    labels.append(f"alpha_{i}->{a}+beta_{j}->{b}")

    return dets, labels


def determinant_overlap(det_bra: Det, det_ket: Det, S_MO: np.ndarray) -> float:
    """Overlap of two Slater determinants given the MO overlap matrix S_MO.

    <det_bra|det_ket> = det(S_MO[alpha_bra, alpha_ket]) * det(S_MO[beta_bra, beta_ket])

    Returns 0 if the alpha (or beta) determinant sizes differ.
    """
    a_bra, b_bra = det_bra
    a_ket, b_ket = det_ket
    if len(a_bra) != len(a_ket) or len(b_bra) != len(b_ket):
        return 0.0

    if len(a_bra) == 0:
        da = 1.0
    else:
        sub_a = S_MO[np.ix_(list(a_bra), list(a_ket))]
        da = np.linalg.det(sub_a)

    if len(b_bra) == 0:
        db = 1.0
    else:
        sub_b = S_MO[np.ix_(list(b_bra), list(b_ket))]
        db = np.linalg.det(sub_b)

    return float(da * db)


def cisd_overlap(ci_bra, dets_bra, ci_ket, dets_ket, S_MO) -> float:
    """Full CISD overlap: <Psi_bra|Psi_ket> = sum_{AB} c_A^bra c_B^ket <det_A|det_B>."""
    total = 0.0
    for c_b, d_b in zip(ci_bra, dets_bra):
        for c_k, d_k in zip(ci_ket, dets_ket):
            total += c_b * c_k * determinant_overlap(d_b, d_k, S_MO)
    return total


def normalize_ci(ci):
    """Normalize a CI vector to unit 2-norm."""
    n = np.linalg.norm(ci)
    return ci / n if n > 0 else ci
