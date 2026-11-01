"""Tests for phase correction."""
import numpy as np
import pytest
from cis_nto.phase import correct_phase_matrix, verify_corrected


def test_identity_no_change():
    """Diagonal identity matrix should not be changed."""
    overlap = np.eye(3, dtype=complex)[None, :, :]
    corrected = correct_phase_matrix(overlap)
    np.testing.assert_allclose(corrected, overlap)


def test_flip_negative_diagonal():
    """Negative diagonal should be flipped to positive."""
    overlap = np.array([[
        [-1.0, 0.1, 0.0],
        [0.1,  1.0, 0.0],
        [0.0,  0.0, 1.0],
    ]], dtype=complex)
    corrected, swaps, flips = correct_phase_matrix(overlap, return_swaps=True)
    assert np.real(corrected[0, 0, 0]) > 0
    assert flips[0] == 1


def test_hungarian_swaps_states():
    """Swapped states 1 and 2 should be reordered."""
    overlap = np.array([[
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],   # state 1 at position 2
        [0.0, 1.0, 0.0],   # state 2 at position 1
    ]], dtype=complex)
    corrected, swaps, flips = correct_phase_matrix(overlap, return_swaps=True)
    assert swaps[0] == 2  # states 1 and 2 swapped
    assert np.allclose(np.diag(corrected[0]), [1.0, 1.0, 1.0])


def test_verify_corrected_passes():
    """Positive diagonal should pass verification."""
    overlap = np.eye(3, dtype=complex)[None, :, :]
    ok, n_bad = verify_corrected(overlap)
    assert ok
    assert n_bad == 0


def test_verify_corrected_fails():
    """Negative diagonal should fail verification."""
    overlap = -np.eye(3, dtype=complex)[None, :, :]
    ok, n_bad = verify_corrected(overlap)
    assert not ok
    assert n_bad == 3


def test_wrong_shape_raises():
    """Should raise on 2D input."""
    with pytest.raises(ValueError):
        correct_phase_matrix(np.eye(3, dtype=complex))
