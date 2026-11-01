"""Verify L4M formula against Phase 1 brute-force."""
import numpy as np
import pytest
import sys
sys.path.insert(0, "/home/hamid/cis_nto/test_cisd")
from bruteforce_cisd import cisd_overlap, normalize_ci, make_determinant
from l4m_cisd import cisd_overlap_l4m


def _setup_3x3():
    S = np.eye(3); S[0, 1] = S[1, 0] = 0.05
    return S


def test_ref_only():
    S = _setup_3x3()
    dets = [make_determinant([0], [0])]
    ci = np.array([1.0])
    brute = cisd_overlap(ci, dets, ci, dets, S)
    l4m = cisd_overlap_l4m(ci, ci, dets, dets, S, S, n_occ_ref=1)
    print(f"ref-only: brute={brute:.15f} l4m={l4m:.15f} diff={abs(brute-l4m):.2e}")
    np.testing.assert_allclose(brute, l4m, atol=1e-13)


def test_single():
    S = _setup_3x3()
    dets = [make_determinant([0], [0]), make_determinant([1], [0])]
    ci = normalize_ci(np.array([0.9, 0.3]))
    brute = cisd_overlap(ci, dets, ci, dets, S)
    l4m = cisd_overlap_l4m(ci, ci, dets, dets, S, S, n_occ_ref=1)
    print(f"single: brute={brute:.15f} l4m={l4m:.15f} diff={abs(brute-l4m):.2e}")
    np.testing.assert_allclose(brute, l4m, atol=1e-13)


def test_full_cisd_3x3():
    S = _setup_3x3()
    dets = [make_determinant([0],[0]), make_determinant([1],[0]),
            make_determinant([0],[1]), make_determinant([1],[1])]
    ci = normalize_ci(np.array([0.9, 0.3, 0.2, 0.1]))
    brute = cisd_overlap(ci, dets, ci, dets, S)
    l4m = cisd_overlap_l4m(ci, ci, dets, dets, S, S, n_occ_ref=1)
    print(f"full-CISD: brute={brute:.15f} l4m={l4m:.15f} diff={abs(brute-l4m):.2e}")
    np.testing.assert_allclose(brute, l4m, atol=1e-13)


def test_two_doubles():
    S = _setup_3x3()
    dets = [make_determinant([0],[0]), make_determinant([1],[1])]
    ci = normalize_ci(np.array([0.8, 0.6]))
    brute = cisd_overlap(ci, dets, ci, dets, S)
    l4m = cisd_overlap_l4m(ci, ci, dets, dets, S, S, n_occ_ref=1)
    print(f"ref+double: brute={brute:.15f} l4m={l4m:.15f} diff={abs(brute-l4m):.2e}")
    np.testing.assert_allclose(brute, l4m, atol=1e-13)


def test_larger_6x6():
    n = 6
    S = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            S[i,j] = S[j,i] = 0.02
    dets = [make_determinant([0,1],[0,1]),
            make_determinant([2,3],[0,1]),
            make_determinant([0,1],[2,3]),
            make_determinant([2,3],[2,3])]
    ci = normalize_ci(np.array([0.85, 0.25, 0.25, 0.10]))
    brute = cisd_overlap(ci, dets, ci, dets, S)
    l4m = cisd_overlap_l4m(ci, ci, dets, dets, S, S, n_occ_ref=2)
    print(f"6x6: brute={brute:.15f} l4m={l4m:.15f} diff={abs(brute-l4m):.2e}")
    np.testing.assert_allclose(brute, l4m, atol=1e-13)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
