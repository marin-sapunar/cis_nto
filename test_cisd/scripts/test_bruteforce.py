"""Tests for brute-force CISD overlap."""
import numpy as np
import sys
sys.path.insert(0, "/home/hamid/cis_nto/test_cisd")
from bruteforce_cisd import (
    make_determinant, gen_cisd_determinants,
    determinant_overlap, cisd_overlap, normalize_ci,
)


def test_self_overlap_normalized():
    """Normalized CISD wavefunction self-overlap should be exactly 1."""
    # 2 occupied, 1 virtual: gives 1 ref + 2 singles + 1 double = 4 dets
    dets, labels = gen_cisd_determinants(n_occ=1, n_virt=1)
    # Actually: n_occ=1, n_virt=1 -> ref + 1 alpha single + 1 beta single + 1 double = 4
    S_MO = np.eye(2)  # identity overlap
    ci = np.array([0.8, 0.3, 0.4, 0.2])
    ci = normalize_ci(ci)
    S = cisd_overlap(ci, dets, ci, dets, S_MO)
    print(f"Self-overlap: {S:.14f}")
    assert abs(S - 1.0) < 1e-12, f"Expected 1.0, got {S}"


def test_reference_only():
    """A pure reference wavefunction has overlap 1.0 with itself."""
    dets, _ = gen_cisd_determinants(n_occ=1, n_virt=1)
    ci = np.array([1.0, 0.0, 0.0, 0.0])
    S_MO = np.eye(2)
    S = cisd_overlap(ci, dets, ci, dets, S_MO)
    assert abs(S - 1.0) < 1e-12


def test_identity_orthogonal_states():
    """Two different determinants should have ~0 overlap with identity S_MO."""
    # Use alpha single vs beta single: they overlap in the OTHER spin channel
    # ref = (0,),(0,); alpha_single = (1,),(0,); beta_single = (0,),(1,)
    ref = make_determinant([0], [0])
    a_single = make_determinant([1], [0])
    S_MO = np.eye(2)
    # Overlap: alpha channel 0 -> 1 = 0, beta channel 0 -> 0 = 1, so 0*1 = 0
    assert abs(determinant_overlap(a_single, ref, S_MO)) < 1e-12


def test_diagonal_physical():
    """Physical MO overlaps (near-identity) should give self-overlap near 1."""
    # 2 orbitals, small mixing
    S_MO = np.array([
        [1.0,  0.05],
        [0.05, 1.0 ],
    ])
    dets, _ = gen_cisd_determinants(n_occ=1, n_virt=1)
    ci = np.array([0.9, 0.3, 0.2, 0.1])
    ci = normalize_ci(ci)
    S = cisd_overlap(ci, dets, ci, dets, S_MO)
    print(f"Physical S_MO, self-overlap: {S:.10f}")
    assert abs(S - 1.0) < 0.15  # Allow some deviation from MO overlap drift


def test_pure_double_excitation():
    """Test that a pure double excitation has overlap with the reference determined by MO overlaps^2."""
    S_MO = np.array([
        [1.0,  0.1],
        [0.1,  1.0 ],
    ])
    ref = make_determinant([0], [0])
    double = make_determinant([1], [1])
    # <double|ref> = S_MO[1,0] * S_MO[1,0] = 0.1 * 0.1 = 0.01
    ov = determinant_overlap(double, ref, S_MO)
    expected = 0.1 * 0.1
    print(f"<double|ref>: {ov:.10f}, expected: {expected:.10f}")
    assert abs(ov - expected) < 1e-12


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v", "--tb=short"])
